#include"TSinkhornKernel.h"

TSinkhornKernelGenerator::TSinkhornKernelGenerator(
		THierarchicalPartition *_HPX, THierarchicalPartition *_HPY,
		const double * _muX, const double *_muY,
		THierarchicalCostFunctionProvider *_costProvider,
		int _layerBottom
		) {

	HPX=_HPX;
	HPY=_HPY;
	if((_muX!=NULL) && (_muY!=NULL)) {
		muX=_muX;
		muY=_muY;
		useReferenceMeasures=true;
	} else {
		muX=NULL;
		muY=NULL;
		useReferenceMeasures=false;
	}
	costProvider=_costProvider;
	layerBottom=_layerBottom;
	
	useSafeMode=false;
	useFixDuals=false;
	
}


TKernelMatrix TSinkhornKernelGenerator::generate(double _eps, double _slack, bool safeMode, bool fixDuals) {
	
	eps=_eps;
	slack=_slack;

	// dimensions of kernel at given layer
	int xres=HPX->layers[layerBottom]->nCells;
	int yres=HPY->layers[layerBottom]->nCells;
	
	// delete all entries (just in case not yet empty)
	entries.clear();
	
	// check all cells at coarsests layer, which then recursively does hierarchical search
	int nX=HPX->layers[0]->nCells;
	int nY=HPY->layers[0]->nCells;
	for(int iX=0;iX<nX;iX++) {
		for(int iY=0;iY<nY;iY++) {
			checkCell(0,iX,iY);
		}
	}

	// safe mode: detect empty cols and rows and fix them
	std::vector<bool> foundRow;
	std::vector<bool> foundCol;
	if(safeMode) {
		// indicators, which col or row has at least one entry
		foundRow.resize(xres,false);
		foundCol.resize(yres,false);
		// go through all elements and mark corresponding cols and rows
		for(uint i=0;i<entries.size();i++) {
			foundRow[entries[i].row()]=true;
			foundCol[entries[i].col()]=true;
		}
		
		// now go through rows and fix empty rows
		for(int x=0;x<xres;x++) {
			if(!foundRow[x]) {
				eprintf("\tfixing row %d\n",x);
				// find suitable kernel entries
				std::vector<TSinkhornKernelGenerator::TCandidate> data=findKernelLine(x,0);
				
				double offset=0;
				if(fixDuals) {
					offset=data[0].v; // offset value of smallest new entry
					costProvider->alpha[layerBottom][x]+=offset; // adjust corresponding value of alpha
				}
				
				// add new entries to list
				for(uint i=0;i<data.size();i++) {
					foundCol[data[i].z]=true; // mark corresponding col as non-empty
					addEntry(x,data[i].z,data[i].v-offset); // add offset value to entry list
				}
			}
		}
		if(fixDuals) {
			// propagate hierarchical alpha
			HPX->signal_propagate_double(costProvider->alpha, 0, layerBottom, THierarchicalPartition::MODE_MAX);
		}
				
		// now fix empty cols
		for(int y=0;y<yres;y++) {
			if(!foundCol[y]) {
				eprintf("\tfixing col %d\n",y);
				// find suitable kernel entries
				std::vector<TSinkhornKernelGenerator::TCandidate> data=findKernelLine(y,1);
				
				double offset=0;
				if(fixDuals) {
					offset=data[0].v; // offset value of smallest new entry
					costProvider->beta[layerBottom][y]+=offset; // adjust corresponding value of beta
				}
				
				
				// add new entries to list
				for(uint i=0;i<data.size();i++) {
					addEntry(data[i].z,y,data[i].v-offset); // add offset value to entry list
				}
			}
		}
		if(fixDuals) {
			// propagate hierarchical beta
			HPY->signal_propagate_double(costProvider->beta, 0, layerBottom, THierarchicalPartition::MODE_MAX);
		}		
	}

	// assemble kernel matrix
	TKernelMatrix result(xres,yres);
	result.setFromTriplets(entries.begin(), entries.end());
	
	// delete all entries
	entries.clear();	
	return result;
}

void TSinkhornKernelGenerator::checkCell(int layer, int x, int y) {
	// compute lower bound for effective cost function value at given level and entry
	double value=costProvider->getCostEff(layer, x, y);
	if(value<=slack) {
		// value is too small to ignore
		if(layer==layerBottom) {
			// if at finest layer: add entry
			addEntry(x,y,value);
		} else {
			// check finer cells
			int nChildrenX=HPX->layers[layer]->nChildren[x];
			int nChildrenY=HPY->layers[layer]->nChildren[y];
			int *childrenX=HPX->layers[layer]->children[x];
			int *childrenY=HPY->layers[layer]->children[y];
			for(int iX=0;iX<nChildrenX;iX++) {
				for(int iY=0;iY<nChildrenY;iY++) {
					checkCell(layer+1,childrenX[iX],childrenY[iY]);
				}
			}
		}
	}
	
}


TKernelMatrix TSinkhornKernelGenerator::refine(const TKernelMatrix &oldKernel, double _eps) {
	// it is assumed that oldKernel is a kernel matrix on layerBottom-1
	eps=_eps;
	entries.clear();
	
	// iterate over elements of old kernel
	for (int xOld=0; xOld<oldKernel.outerSize(); xOld++) {
		for (TKernelMatrix::InnerIterator it(oldKernel,xOld); it; ++it) {
			int yOld=it.col();
			// for each element of old kernel compute refined elements
			// check finer cells
			int nChildrenX=HPX->layers[layerBottom-1]->nChildren[xOld];
			int nChildrenY=HPY->layers[layerBottom-1]->nChildren[yOld];
			int *childrenX=HPX->layers[layerBottom-1]->children[xOld];
			int *childrenY=HPY->layers[layerBottom-1]->children[yOld];
			for(int iX=0;iX<nChildrenX;iX++) {
				for(int iY=0;iY<nChildrenY;iY++) {
					addEntry(childrenX[iX],childrenY[iY]);
				}
			}
		}
	}
		
	// assemble kernel matrix
	TKernelMatrix result(HPX->layers[layerBottom]->nCells,HPY->layers[layerBottom]->nCells);
	result.setFromTriplets(entries.begin(), entries.end());
	entries.clear();
	return result;
}



std::vector<TSinkhornKernelGenerator::TCandidate> TSinkhornKernelGenerator::findKernelLine(int a, int mode) {
	// mode=0: fix row with x=a
	// else: fix col with y=a

	std::vector<TSinkhornKernelGenerator::TCandidate> result(0);
	bool foundSmallestElement=false; // indicates whether smallest element in kernel row has already
			// been found in hierarchical search
	double smallestValue=0; // this stores the corresponding value, once smallest element was found

	// deciding whether to fix row or col:
	// set HPB to secondary partition (e.g. HPY if row mode)
	THierarchicalPartition *HPB;
	if(mode==0) {
		HPB=HPY;
	} else {
		HPB=HPX;
	}
	// for queries to costProvider use simple if
	
	// generate hierarchical kernel row at coarsest level
	int nB=HPB->layers[0]->nCells;
	std::vector<TSinkhornKernelGenerator::TCandidate> data(nB);
	for(int iB=0;iB<nB;iB++) {
		data[iB].layer=0;
		data[iB].z=iB;
		if(mode==0) {
			data[iB].v=costProvider->getCostEffAsym(layerBottom, a, 0, iB); // new asymmetric call, where x is always at finest level
		} else {
			data[iB].v=costProvider->getCostEffAsym(0, iB, layerBottom, a); // new asymmetric call, where x is always at finest level
		}
	}
	
	// create initial kernel row
	TSinkhornKernelGenerator::TCandidateList row(data.data(),nB);
	
	// now run recursive refinement strategy
	while((row.data.begin() != row.data.end())) {
		// extract smallest element from front of list
		TSinkhornKernelGenerator::TCandidate smallestEntry=row.data.front();
		// afterwards, remove element from list
		row.data.pop_front();

		//eprintf("%d\t%d\t%e\n",smallestEntry.layer,smallestEntry.z,smallestEntry.v);

		if(smallestEntry.layer<layerBottom) {
			// refinement required
			int newLayer=smallestEntry.layer+1; // new layer
			int nChildrenB=HPB->layers[smallestEntry.layer]->nChildren[smallestEntry.z]; // number of children of b on new layer
			int *childrenB=HPB->layers[smallestEntry.layer]->children[smallestEntry.z]; // pointer to b children

			// for children of b obtain new entries for hierarchical kernel row
			data.resize(nChildrenB);			
			for(int iB=0;iB<nChildrenB;iB++) {
				data[iB].layer=newLayer;
				data[iB].z=childrenB[iB];
				if(mode==0) {
					data[iB].v=costProvider->getCostEffAsym(layerBottom, a, newLayer, childrenB[iB]); // see above
				} else {
					data[iB].v=costProvider->getCostEffAsym(newLayer, childrenB[iB], layerBottom, a); // see above
				}
				//eprintf("\t%d\t%d\t%e\n",data[iY].layer,data[iY].z,data[iY].v);
			}
			
			// merge these new entries into list at appropriate locations
			row.merge(data.data(),nChildrenB);			
			
		} else {
			// front element of hiearchical kernel row is now on finest level

			if(!foundSmallestElement) {
				// if found first element
				// add as first element to result list, set "smallestValue" accordingly
				foundSmallestElement=true;
				smallestValue=smallestEntry.v;
				result.push_back(smallestEntry);
			} else {
				// if already found smallest element
				// check slack relative to smallestValue, if sufficiently small, add to result list, else return
				if(smallestEntry.v-smallestValue<=slack) {
					result.push_back(smallestEntry);
				} else {
					return result;
				}
			}
			
		}
	}
	// should the finest element never exceed the slack value, simply return full list now
	return result;	
	
}





//THierarchicalKernelRow::THierarchicalKernelRow(THierarchicalKernelRowEntry *newData, int n) {
//	// construct linked list from array data
//	
//	// first sort
//	std::sort(newData,newData+n);
//	
//	// then iteratively add elements
//	std::forward_list<THierarchicalKernelRowEntry>::iterator it=data.before_begin();
//	for(int i=0;i<n;i++) {
//		it=data.insert_after(it,newData[i]);
//	}
//	// set number of elements
//	size=n;
//}

//THierarchicalKernelRow::~THierarchicalKernelRow() {
//}

//void THierarchicalKernelRow::merge(THierarchicalKernelRowEntry *newData, int n) {
//	// first sort
//	std::sort(newData,newData+n);

//	// if current list is empty, simply add new elements
//	if(size==0) {
//		// then iteratively add elements
//		std::forward_list<THierarchicalKernelRowEntry>::iterator it=data.before_begin();
//		for(int i=0;i<n;i++) {
//			it=data.insert_after(it,newData[i]);
//		}
//		// set size
//		size=n;
//		return;
//	}
//	
//	// keep two iterators pointing to subsequent elements
//	std::forward_list<THierarchicalKernelRowEntry>::iterator it=data.before_begin();
//	std::forward_list<THierarchicalKernelRowEntry>::iterator itNext=data.begin();

//	int i=0;
//	// go through list, as long as new elements remain to be added
//	while((itNext!=data.end()) && (i<n)) {
//		if(newData[i]<*itNext) {
//			// if new element should be between it and itNext
//			
//			// insert element after it, set it to new element
//			it=data.insert_after(it,newData[i]);
//			// set itNext to subsequent element
//			itNext=it;
//			itNext++;
//			// focus on next new element
//			i++;
//		} else {
//			// otherwise, simply go through list
//			it++;
//			itNext++;
//		}
//	}
//	while(i<n) {
//		// if some new elements are still left, add them to end of list
//		it=data.insert_after(it,newData[i]);
//		i++;
//	}
//	// set new size
//	size+=n;
//	
//}


void SinkhornAbsorbScaling(
		THierarchicalPartition *HP,
		double **alpha,
		TMarginalVector &u,
		int layer,
		double eps
		) {
	
	// absorb values of u into finest layer of alpha
	int xres=HP->layers[layer]->nCells;
	for(int x=0;x<xres;x++) {
		alpha[layer][x]+=eps*log(u[x]);
//		// DEBUG
//		if(!std::isfinite(alpha[layer][x])) {
//			eprintf("NaN at pos %d\t%f\t%f\n",x,alpha[layer][x],u[x]);
//		}
//		// END DEBUG
		u[x]=1.;
	}
	// propagate maximal values over hierarchical alpha values
	HP->signal_propagate_double(alpha, 0, layer, THierarchicalPartition::MODE_MAX);
	
}

TSparseCSRContainer SinkhornKernelGetCSRData(const TKernelMatrix& kernel) {
	TSparseCSRContainer result;
	result.xres=kernel.outerSize();
	result.yres=kernel.innerSize();
	result.nonZeros=kernel.nonZeros();
	result.data.resize(result.nonZeros);
	result.indices.resize(result.nonZeros);
	result.indptr.resize(result.xres+1);
	

	// iterate over elements of kernel
	int index=0; // counter over all indices to set data,indices and indptr appropriately
	result.indptr[0]=0;
	for (int x=0; x<kernel.outerSize(); x++) {
		for (TKernelMatrix::InnerIterator it(kernel,x); it; ++it) {
			result.data[index]=it.value();
			result.indices[index]=it.col();
			index++;
		}
		result.indptr[x+1]=index;
	}
	
	return result;


}

TSparsePosContainer SinkhornKernelGetPosData(const TKernelMatrix& kernel) {

	TSparsePosContainer result;
	result.xres=kernel.outerSize();
	result.yres=kernel.innerSize();
	result.nParticles=kernel.nonZeros();
	result.mass.resize(result.nParticles);
	result.posStart.resize(result.nParticles);
	result.posEnd.resize(result.nParticles);
	

	// iterate over elements of kernel
	int index=0; // counter over all indices
	for (int x=0; x<kernel.outerSize(); x++) {
		for (TKernelMatrix::InnerIterator it(kernel,x); it; ++it) {
			result.mass[index]=it.value();
			result.posStart[index]=it.row();
			result.posEnd[index]=it.col();
			index++;
		}
	}
	
	return result;

}


template class THierarchicalSearchList<TSinkhornKernelGenerator::TCandidate>;

