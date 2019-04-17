#include"MultiScaleTools.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TMultiScaleSetupBase
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//TMultiScaleSetupBase::TMultiScaleSetupBase(TDoubleMatrix *_posX, TDoubleMatrix *_posY, double *_muX, double *_muY,
//		int _depth) {
//	posX=_posX;
//	posY=_posY;
//	muX=_muX;
//	muY=_muY;
//	depth=_depth;
//	
//	HBX=NULL;
//	HBY=NULL;
//	
//	alphaH=NULL;
//	betaH=NULL;
//	xRadii=NULL;
//	yRadii=NULL;

//	xNeighboursH=NULL;

//	HierarchyBuilderChildMode=THierarchyBuilder::CM_Grid;

//}

//TMultiScaleSetupBase::~TMultiScaleSetupBase() {
//	// free dual variables if allocated
//	if(alphaH!=NULL) {
//		HPX->signal_free_double(alphaH, 0, HPX->nLayers-1);
//	}
//	if(betaH!=NULL) {
//		HPY->signal_free_double(betaH, 0, HPY->nLayers-1);
//	}

//	// free radii if allocated
//	if(xRadii!=NULL) {
//		HPX->signal_free_double(xRadii, 0, HPX->nLayers-2);
//	}
//	if(yRadii!=NULL) {
//		HPY->signal_free_double(yRadii, 0, HPY->nLayers-2);
//	}

//	// free neighbours if assigned
//	if(xNeighboursH!=NULL) {
//		for(int i=0;i<nLayers;i++) {
//			delete xNeighboursH[i];
//		}
//		free(xNeighboursH);
//	}


//	// free hiearchical partitions

//	if(HBX!=NULL) {

//		free(xresH);		
//		HBX->freeSignal(muXH,HBX->layers.size());
//		HBX->freeSignal(posXH,HBX->layers.size());
//		delete HPX;
//		delete HBX;
//	}
//	if(HBY!=NULL) {

//		free(yresH);		
//		HBY->freeSignal(muYH,HBY->layers.size());
//		HBY->freeSignal(posYH,HBY->layers.size());
//		delete HPY;
//		delete HBY;
//	
//		
//	}	
//	

//}

//int TMultiScaleSetupBase::BasicSetup() {
//	if((posX->depth!=2) || (posY->depth!=2)) {
//		eprintf("ERROR: marginal point clouds must be 2d arrays.\n");
//		return ERR_PREP_INIT_POINTCLOUD2D;
//	}
//	if(posX->dimensions[1]!=posY->dimensions[1]) {
//		return ERR_PREP_INIT_DIMMISMATCH;
//	}
//	xres=posX->dimensions[0];
//	yres=posY->dimensions[0];
//	dim=posX->dimensions[1];
//	return 0;
//}

//int TMultiScaleSetupBase::BasicMeasureChecks() {

//	// sanity check: muX and muY must be strictly positive
//	if(doubleArrayMin(muX,xres)<=0.) {
//		eprintf("ERROR: minimum of muX is not strictly positive.\n");
//		return ERR_PREP_INIT_MUXNEG;
//	}
//	if(doubleArrayMin(muY,yres)<=0.) {
//		eprintf("ERROR: minimum of muY is not strictly positive.\n");
//		return ERR_PREP_INIT_MUYNEG;
//	}

//	// TODO: more sanity checks?

//	return 0;
//}


//int TMultiScaleSetupBase::SetupHierarchicalPartition(double *mu, double *pos, int res, int dim, int depth,
//		THierarchyBuilder **_HB, THierarchicalPartition **_HP, double ***_posH, double ***_muH, int **_resH) {
//	
//	// create hierarchical partition
//	THierarchyBuilder *HB = new THierarchyBuilder(pos,res,dim, HierarchyBuilderChildMode, depth);
//		
//	// convert to format used by hierarchical solver
//	THierarchicalPartition *HP=HB->convert();
//	
//	// hierarchical positions
//	// position of nodes in coarser hierarchy levels
//	double **posH=HB->allocateDoubleSignal(dim);
//	HB->getSignalPos(posH);
//	
//	// hierarchical masses
//	// combined masses of nodes in coarser hierarchy levels
//	double **muH=HB->allocateDoubleSignal(1);
//	HP->computeHierarchicalMasses(mu,muH);

//	// get hierarchical cardinalities of each marginal on each level
//	int *resH=HB->getResH();

//	*_HB=HB;
//	*_HP=HP;
//	*_posH=posH;
//	*_muH=muH;
//	*_resH=resH;
//	
//	return 0;
//	
//}

//int TMultiScaleSetupBase::SetupHierarchicalPartitions() {
//	int msg;
//	msg=SetupHierarchicalPartition(muX, posX->data, xres, dim, depth,
//			&HBX, &HPX, &posXH, &muXH, &xresH);
//	if(msg!=0) {
//		return msg;
//	}
//	
//	msg=SetupHierarchicalPartition(muY, posY->data, yres, dim, depth,
//			&HBY, &HPY, &posYH, &muYH, &yresH);
//	if(msg!=0) {
//		return msg;
//	}
//			
//	nLayers=HPX->nLayers;

//	return 0;
//}

//int TMultiScaleSetupBase::Setup() {
//	int msg;

//	msg=BasicSetup();
//	if(msg!=0) { return msg; }
//	msg=BasicMeasureChecks();
//	if(msg!=0) { return msg; }
//	msg=SetupHierarchicalPartitions();
//	if(msg!=0) { return msg; }	

//	// only invoke these on demand
////	msg=SetupDuals();
////	if(msg!=0) { return msg; }	
////	msg=SetupRadii();
////	if(msg!=0) { return msg; }	


//	return 0;

//}


//int TMultiScaleSetupBase::SetupDuals() {
//	alphaH=HPX->signal_allocate_double(0,HPX->nLayers-1);
//	betaH=HPY->signal_allocate_double(0,HPY->nLayers-1);
//	return 0;
//}

//int TMultiScaleSetupBase::SetupRadii() {
//	xRadii=HBX->getSignalRadii();	
//	yRadii=HBY->getSignalRadii();
//	return 0;
//}


////////////////////////////////////////////
//// Cartesian Grid
////////////////////////////////////////////

//TMultiScaleSetupGrid::TMultiScaleSetupGrid(TDoubleMatrix *_muXGrid, TDoubleMatrix *_muYGrid, int _depth) :
//	TMultiScaleSetupBase(NULL, NULL, _muXGrid->data, _muYGrid->data, _depth) {
//	muXGrid=_muXGrid;
//	muYGrid=_muYGrid;
//	posX=GridToolsGetGridMatrix(muXGrid->depth,muXGrid->dimensions);
//	posY=GridToolsGetGridMatrix(muYGrid->depth,muYGrid->dimensions);
//	xDimH=NULL;
//	yDimH=NULL;
//		
//}

//int TMultiScaleSetupGrid::SetupHierarchicalPartitions() {
//	int msg;
//	
//	msg=TMultiScaleSetupBase::SetupHierarchicalPartitions(); // call original method
//	if(msg!=0) {
//		return msg;
//	}
//	
//	// initialize xDimH and yDimH
//	xDimH=HBX->getDimH(muXGrid->dimensions);
//	yDimH=HBY->getDimH(muYGrid->dimensions);
//	return 0;
//}

//int TMultiScaleSetupGrid::SetupGridNeighboursX() {
//	if(xNeighboursH!=NULL) {
//		for(int i=0;i<nLayers;i++) {
//			delete xNeighboursH[i];
//		}
//		free(xNeighboursH);
//	}

//	xNeighboursH=(TVarListHandler**) malloc(sizeof(TVarListHandler*)*nLayers);
//	for(int i=0;i<nLayers;i++) {
//		xNeighboursH[i]=new TVarListHandler();
//		xNeighboursH[i]->setupEmpty(xresH[i]);
//		GridToolsGetNeighbours(dim, xDimH+i*dim, xNeighboursH[i]);
//	}
//	
//	return 0;
//}


//TMultiScaleSetupGrid::~TMultiScaleSetupGrid() {
//	freeTDoubleMatrix(posX);
//	freeTDoubleMatrix(posY);
//	
//	if(xDimH!=NULL) {
//		free(xDimH);
//	}
//	if(yDimH!=NULL) {
//		free(yDimH);
//	}
//	
//}


//////////////////////////////////////////
// Barycenter
//////////////////////////////////////////


TMultiScaleSetupBarycenterBase::TMultiScaleSetupBarycenterBase(
		int _nMarginals,
		TDoubleMatrix **_pos, TDoubleMatrix *_posZ,
		double **_mu, double *_muZ,
		int _depth) {

	nMarginals=_nMarginals;
	pos=_pos;
	posZ=_posZ;
	mu=_mu;
	muZ=_muZ;
	depth=_depth;
	
	res=NULL;
	
	HB=NULL;
	HBZ=NULL;
	
	alphaH=NULL;
	betaH=NULL;
	radii=NULL;
	zRadii=NULL;


	HierarchyBuilderChildMode=THierarchyBuilder::CM_Grid;

}

TMultiScaleSetupBarycenterBase::~TMultiScaleSetupBarycenterBase() {
	// free dual variables if allocated
	if(alphaH!=NULL) {
		for(int i=0;i<nMarginals;i++) {
			HP[i]->signal_free_double(alphaH[i], 0, HP[i]->nLayers-1);
		}
		free(alphaH);
	}
	if(betaH!=NULL) {
		for(int i=0;i<nMarginals;i++) {
			HPZ->signal_free_double(betaH[i], 0, HPZ->nLayers-1);
		}
		free(betaH);
	}

	// free radii if allocated
	if(radii!=NULL) {
		for(int i=0;i<nMarginals;i++) {
			HP[i]->signal_free_double(radii[i], 0, HP[i]->nLayers-2);
		}
		free(radii);
	}
	if(zRadii!=NULL) {
		HPZ->signal_free_double(zRadii, 0, HPZ->nLayers-2);
	}

	if(HB!=NULL) {

		for(int i=0;i<nMarginals;i++) {
			free(resH[i]);
			HB[i]->freeSignal(muH[i],HB[i]->layers.size());
			HB[i]->freeSignal(posH[i],HB[i]->layers.size());
			delete HP[i];
			delete HB[i];
		}
		free(muH);
		free(posH);
		free(resH);
		free(HP);
		free(HB);
		
	}
	if(HBZ!=NULL) {

		free(zresH);		
		HBZ->freeSignal(muZH,HBZ->layers.size());
		HBZ->freeSignal(posZH,HBZ->layers.size());
		delete HPZ;
		delete HBZ;
	
		
	}
	
	if(res!=NULL) {
		free(res);
	}
	
}


int TMultiScaleSetupBarycenterBase::BasicSetup() {

	res=(int*) malloc(sizeof(int)*nMarginals);

	for(int i=0;i<nMarginals;i++) {
		if(pos[i]->depth!=2) {
			eprintf("ERROR: marginal point clouds must be 2d arrays.\n");
			return ERR_PREP_INIT_POINTCLOUD2D;
		}
		if(pos[i]->dimensions[1]!=posZ->dimensions[1]) {
			return ERR_PREP_INIT_DIMMISMATCH;
		}
		res[i]=pos[i]->dimensions[0];
	}
	if(posZ->depth!=2) {
		eprintf("ERROR: marginal point clouds must be 2d arrays.\n");
		return ERR_PREP_INIT_POINTCLOUD2D;
	}
	zres=posZ->dimensions[0];
	dim=posZ->dimensions[1];
	return 0;
}

int TMultiScaleSetupBarycenterBase::BasicMeasureChecks() {

	// sanity check: measures must be strictly positive
	for(int i=0;i<nMarginals;i++) {
		if(doubleArrayMin(mu[i],res[i])<=0.) {
			eprintf("ERROR: minimum of muX is not strictly positive.\n");
			return ERR_PREP_INIT_MUXNEG;
		}
	}
	if(doubleArrayMin(muZ,zres)<=0.) {
		eprintf("ERROR: minimum of muY is not strictly positive.\n");
		return ERR_PREP_INIT_MUYNEG;
	}

	// TODO: more sanity checks?

	return 0;
}

int TMultiScaleSetupBarycenterBase::SetupHierarchicalPartition(double *aMu, double *aPos, int aRes, int dim, int depth,
		THierarchyBuilder **_aHB, THierarchicalPartition **_aHP, double ***_aPosH, double ***_aMuH, int **_aResH) {
	
	// create hierarchical partition
	THierarchyBuilder *aHB = new THierarchyBuilder(aPos,aRes,dim, HierarchyBuilderChildMode, depth);
		
	// convert to format used by hierarchical solver
	THierarchicalPartition *aHP=aHB->convert();
	
	// hierarchical positions
	// position of nodes in coarser hierarchy levels
	double **aPosH=aHB->allocateDoubleSignal(dim);
	aHB->getSignalPos(aPosH);
	
	// hierarchical masses
	// combined masses of nodes in coarser hierarchy levels
	double **aMuH=aHB->allocateDoubleSignal(1);
	aHP->computeHierarchicalMasses(aMu,aMuH);

	// get hierarchical cardinalities of each marginal on each level
	int *aResH=aHB->getResH();

	*_aHB=aHB;
	*_aHP=aHP;
	*_aPosH=aPosH;
	*_aMuH=aMuH;
	*_aResH=aResH;
	
	return 0;
	
}

int TMultiScaleSetupBarycenterBase::SetupHierarchicalPartitions() {
	int msg;
	HB=(THierarchyBuilder**) malloc(sizeof(THierarchyBuilder*)*nMarginals);
	HP=(THierarchicalPartition**) malloc(sizeof(THierarchicalPartition*)*nMarginals);
	posH=(double***) malloc(sizeof(double**)*nMarginals);
	muH=(double***) malloc(sizeof(double**)*nMarginals);
	resH=(int**) malloc(sizeof(int*)*nMarginals);
	
	for(int i=0;i<nMarginals;i++) {
		msg=SetupHierarchicalPartition(mu[i], pos[i]->data, res[i], dim, depth,
				HB+i, HP+i, posH+i, muH+i, resH+i);
		if(msg!=0) { return msg; }
	}
		
	msg=SetupHierarchicalPartition(muZ, posZ->data, zres, dim, depth,
			&HBZ, &HPZ, &posZH, &muZH, &zresH);
	if(msg!=0) { return msg; }

	nLayers=HPZ->nLayers;
		
	return 0;
}

int TMultiScaleSetupBarycenterBase::Setup() {
	int msg;

	msg=BasicSetup();
	if(msg!=0) { return msg; }
	msg=BasicMeasureChecks();
	if(msg!=0) { return msg; }
	msg=SetupHierarchicalPartitions();
	if(msg!=0) { return msg; }	

//	msg=SetupDuals();
//	if(msg!=0) { return msg; }	
//	msg=SetupRadii();
//	if(msg!=0) { return msg; }	


	return 0;

}


int TMultiScaleSetupBarycenterBase::SetupDuals() {
	alphaH=(double***) malloc(sizeof(double**)*nMarginals);
	betaH=(double***) malloc(sizeof(double**)*nMarginals);
	
	for(int i=0;i<nMarginals;i++) {
		alphaH[i]=HP[i]->signal_allocate_double(0,HP[i]->nLayers-1);
		betaH[i]=HPZ->signal_allocate_double(0,HPZ->nLayers-1);
	}
	return 0;
}

int TMultiScaleSetupBarycenterBase::SetupRadii() {
	radii=(double***) malloc(sizeof(double**)*nMarginals);
	for(int i=0;i<nMarginals;i++) {
		radii[i]=HB[i]->getSignalRadii();	
	}
	
	zRadii=HBZ->getSignalRadii();
	return 0;
}


//////////////////////////////////////////
// NN Search
//////////////////////////////////////////


std::vector<int> THierarchicalNN::find(double *posX, double **posYH, double **radiiY,
		THierarchicalPartition *HPY, int layerBottom, int nElements) {

	int dim=HPY->dim;

	std::vector<int> result(0);

	int foundElements=0;
	// initialize search list with only one element at coarsest level
	TCandidate candidate;
	candidate.layer=0;
	candidate.z=0;
	candidate.dist=0.;
	
	TCandidateList candidateList(&candidate,1);

	// list where to collect refinement lists
	std::vector<TCandidate> data(0);
	
	// now run recursive refinement strategy
	while((candidateList.data.begin() != candidateList.data.end())) {
		// extract smallest element from front of list
		TCandidate smallestEntry=candidateList.data.front();
		// afterwards, remove element from list
		candidateList.data.pop_front();

		//eprintf("%d\t%d\t%e\n",smallestEntry.layer,smallestEntry.z,smallestEntry.dist);

		if(smallestEntry.layer<layerBottom) {
			//eprintf("refining\n");
			// refinement required
			int newLayer=smallestEntry.layer+1; // new layer

			// number of children of candidate on new layer
			int nChildren=HPY->layers[smallestEntry.layer]->nChildren[smallestEntry.z];
			// pointer to children
			int *children=HPY->layers[smallestEntry.layer]->children[smallestEntry.z];

			// for children obtain new entries
			data.resize(nChildren);			
			for(int i=0;i<nChildren;i++) {
				data[i].layer=newLayer;
				data[i].z=children[i];
				data[i].dist=std::pow(EUCL_lincombSqr(posX, posYH[newLayer]+dim*data[i].z, 1., -1., dim),0.5);
				if(newLayer<layerBottom) {
					data[i].dist-=radiiY[newLayer][data[i].z];
				}
			}
			
			// merge these new entries into list at appropriate locations
			candidateList.merge(data.data(),nChildren);			
			
		} else {
			//eprintf("found finest element\n");
			// front element of candidate list is now on finest level
			// therefore, add to result
			result.push_back(smallestEntry.z);
			foundElements++;
			// return list, when sufficient number of elements is found
			if(foundElements>=nElements) {
				return result;
			}
			
		}
	}
	// should the number of element never exceed nElements, simply return full list now
	return result;
}

TVarListHandler* THierarchicalNN::getNeighbours(double **posXH, double **radiiX,
		THierarchicalPartition *HPX, int layerBottom, int nElements) {

	int xres=HPX->layers[layerBottom]->nCells;
	int dim=HPX->dim;
	
	TVarListHandler *result=new TVarListHandler;
	result->setupEmpty(xres);
	
	for(int x=0;x<xres;x++) {
		// find nElements+1 nearest neighbours of point
		// +1 because nearest neighbour will always be point itself
		std::vector<int> line=find(posXH[layerBottom]+dim*x, posXH, radiiX,
				HPX, layerBottom, nElements+1);
		
		// add elements from 2nd onwards to result
		result->addToLine(x,line.data()+1,line.size()-1);
	}
	
	return result;
}

TVarListHandler** THierarchicalNN::getNeighboursH(double **posXH, double **radiiX,
		THierarchicalPartition *HPX, int nElements) {
	
	int nLayers=HPX->nLayers;
	TVarListHandler **result=(TVarListHandler**) malloc(sizeof(TVarListHandler*)*nLayers);
	for(int layer=0;layer<nLayers;layer++) {
		result[layer]=getNeighbours(posXH, radiiX,
				HPX, layer, nElements);
	}
	return result;
}




void THierarchicalDualMaximizer::getMaxDual(THierarchicalPartition *partitionX, THierarchicalPartition *partitionY,
		double **alpha, double **beta, int layerFine,
		THierarchicalCostFunctionProvider *costProvider,
		int mode) {

	// for one given dual variable, compute the maximal other one, using hierarchical search

	// mode=MODE_ALPHA: compute alpha for given beta
	// else: compute beta for given alpha

	THierarchicalPartition *HPA,*HPB;

	// deciding whether to fix row or col:
	// set HPA to primary partition
	// set HPB to secondary partition (e.g. HPY if row mode)
	if (mode==MODE_ALPHA) {
		HPA=partitionX;
		HPB=partitionY;
	} else {
		HPA=partitionY;
		HPB=partitionX;
	}
	// for queries to costProvider use simple if

	int ares; // number of entries in row or column that need to be computed
	ares=HPA->layers[layerFine]->nCells;

	for (int a=0; a<ares; a++) {

	
		// generate hierarchical row at coarsest level
		int nB=HPB->layers[0]->nCells;
		std::vector<THierarchicalDualMaximizer::TCandidate> data(nB);
		for(int iB=0;iB<nB;iB++) {
			data[iB].layer=0;
			data[iB].z=iB;
			if(mode==MODE_ALPHA) {
				data[iB].v=costProvider->getCostAsym(layerFine, a, 0, iB)-beta[0][iB];
			} else {
				data[iB].v=costProvider->getCostAsym(0, iB, layerFine, a)-alpha[0][iB];
			}
		}
	
		// create initial kernel row
		THierarchicalDualMaximizer::TCandidateList row(data.data(),nB);
	
		// now run recursive refinement strategy
		while((row.data.begin() != row.data.end())) {
			// extract smallest element from front of list
			THierarchicalDualMaximizer::TCandidate smallestEntry=row.data.front();
			// afterwards, remove element from list
			row.data.pop_front();

			if(smallestEntry.layer<layerFine) {
				// refinement required
				int newLayer=smallestEntry.layer+1; // new layer
				int nChildrenB=HPB->layers[smallestEntry.layer]->nChildren[smallestEntry.z]; // number of children of b on new layer
				int *childrenB=HPB->layers[smallestEntry.layer]->children[smallestEntry.z]; // pointer to b children

				// for children of b obtain new entries for hierarchical row
				data.resize(nChildrenB);			
				for(int iB=0;iB<nChildrenB;iB++) {
					data[iB].layer=newLayer;
					data[iB].z=childrenB[iB];
					if(mode==MODE_ALPHA) {
						data[iB].v=costProvider->getCostAsym(layerFine, a, newLayer, childrenB[iB])-beta[newLayer][childrenB[iB]];
					} else {
						data[iB].v=costProvider->getCostAsym(newLayer, childrenB[iB], layerFine, a)-alpha[newLayer][childrenB[iB]];
					}

				}
			
				// merge these new entries into list at appropriate locations
				row.merge(data.data(),nChildrenB);			
			
			} else {
				// front element of hiearchical row is now on finest level
				if(mode==MODE_ALPHA) {
					alpha[layerFine][a]=smallestEntry.v;
				} else {
					beta[layerFine][a]=smallestEntry.v;
				}
				break;
							
			}
		}
	}	

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TMULTISCALESETUPBASE SINGLE MARGINAL

TMultiScaleSetupSingleBase::TMultiScaleSetupSingleBase(TDoubleMatrix *_pos, double *_mu, int _depth) {
	pos=_pos;
	mu=_mu;
	depth=_depth;
	
	HB=NULL;
	
	alphaH=NULL;
	radii=NULL;

	neighboursH=NULL;

	HierarchyBuilderChildMode=THierarchyBuilder::CM_Grid;

}

TMultiScaleSetupSingleBase::~TMultiScaleSetupSingleBase() {
	// free dual variables if allocated
	if(alphaH!=NULL) {
		HP->signal_free_double(alphaH, 0, HP->nLayers-1);
	}

	// free radii if allocated
	if(radii!=NULL) {
		HP->signal_free_double(radii, 0, HP->nLayers-2);
	}

	// free neighbours if assigned
	if(neighboursH!=NULL) {
		for(int i=0;i<nLayers;i++) {
			delete neighboursH[i];
		}
		free(neighboursH);
	}


	// free hiearchical partitions

	if(HB!=NULL) {

		free(resH);		
		HB->freeSignal(muH,HB->layers.size());
		HB->freeSignal(posH,HB->layers.size());
		delete HP;
		delete HB;
	}

}

int TMultiScaleSetupSingleBase::BasicSetup() {
	if(pos->depth!=2) {
		eprintf("ERROR: marginal point clouds must be 2d arrays.\n");
		return ERR_PREP_INIT_POINTCLOUD2D;
	}
	res=pos->dimensions[0];
	dim=pos->dimensions[1];
	return 0;
}

int TMultiScaleSetupSingleBase::BasicMeasureChecks() {

	// sanity check: mu must be strictly positive
	if(doubleArrayMin(mu,res)<=0.) {
		eprintf("ERROR: minimum of mu is not strictly positive.\n");
		return ERR_PREP_INIT_MUXNEG;
	}

	// TODO: more sanity checks?

	return 0;
}


int TMultiScaleSetupSingleBase::SetupHierarchicalPartition() {
	
	// create hierarchical partition
	HB = new THierarchyBuilder(pos->data,res,dim, HierarchyBuilderChildMode, depth);
		
	// convert to format used by hierarchical solver
	HP=HB->convert();
	
	// hierarchical positions
	// position of nodes in coarser hierarchy levels
	posH=HB->allocateDoubleSignal(dim);
	HB->getSignalPos(posH);
	
	// hierarchical masses
	// combined masses of nodes in coarser hierarchy levels
	muH=HB->allocateDoubleSignal(1);
	HP->computeHierarchicalMasses(mu,muH);

	// get hierarchical cardinalities of each marginal on each level
	resH=HB->getResH();

	nLayers=HP->nLayers;
	
	return 0;
	
}

int TMultiScaleSetupSingleBase::Setup() {
	int msg;

	msg=BasicSetup();
	if(msg!=0) { return msg; }
	msg=BasicMeasureChecks();
	if(msg!=0) { return msg; }
	msg=SetupHierarchicalPartition();
	if(msg!=0) { return msg; }	

	// only invoke these on demand
//	msg=SetupDuals();
//	if(msg!=0) { return msg; }	
//	msg=SetupRadii();
//	if(msg!=0) { return msg; }	


	return 0;

}


int TMultiScaleSetupSingleBase::SetupDuals() {
	alphaH=HP->signal_allocate_double(0,HP->nLayers-1);
	return 0;
}

int TMultiScaleSetupSingleBase::SetupRadii() {
	radii=HB->getSignalRadii();	
	return 0;
}


//////////////////////////////////////////
// Cartesian Grid
//////////////////////////////////////////

TMultiScaleSetupSingleGrid::TMultiScaleSetupSingleGrid(TDoubleMatrix *_muGrid, int _depth) :
	TMultiScaleSetupSingleBase(NULL, _muGrid->data, _depth) {
	muGrid=_muGrid;
	pos=GridToolsGetGridMatrix(muGrid->depth,muGrid->dimensions);
	dimH=NULL;
		
}

int TMultiScaleSetupSingleGrid::SetupHierarchicalPartition() {
	int msg;
	
	msg=TMultiScaleSetupSingleBase::SetupHierarchicalPartition(); // call original method
	if(msg!=0) {
		return msg;
	}
	
	// initialize dimH
	dimH=HB->getDimH(muGrid->dimensions);
	return 0;
}

int TMultiScaleSetupSingleGrid::SetupGridNeighbours() {
	if(neighboursH!=NULL) {
		for(int i=0;i<nLayers;i++) {
			delete neighboursH[i];
		}
		free(neighboursH);
	}

	neighboursH=(TVarListHandler**) malloc(sizeof(TVarListHandler*)*nLayers);
	for(int i=0;i<nLayers;i++) {
		neighboursH[i]=new TVarListHandler();
		neighboursH[i]->setupEmpty(resH[i]);
		GridToolsGetNeighbours(dim, dimH+i*dim, neighboursH[i]);
	}
	
	return 0;
}


TMultiScaleSetupSingleGrid::~TMultiScaleSetupSingleGrid() {
	freeTDoubleMatrix(pos);
	
	if(dimH!=NULL) {
		free(dimH);
	}
	
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


TMultiScaleSetupBarycenterContainer::TMultiScaleSetupBarycenterContainer() {
	nMarginals=0;
	HP=NULL;
	HPZ=NULL;
	muH=NULL;
	muZH=NULL;
	alphaH=NULL;
	betaH=NULL;
	costProvider=NULL;
	weights=NULL;
	resH=NULL;
	resZH=NULL;
}

TMultiScaleSetupBarycenterContainer::TMultiScaleSetupBarycenterContainer(const int _nMarginals)
	: TMultiScaleSetupBarycenterContainer() {
	setupEmpty(_nMarginals);
}

TMultiScaleSetupBarycenterContainer::~TMultiScaleSetupBarycenterContainer() {
	cleanup();
}

void TMultiScaleSetupBarycenterContainer::setupEmpty(const int _nMarginals) {
	cleanup();
	nMarginals=_nMarginals;
	HP=(THierarchicalPartition**) malloc(sizeof(THierarchicalPartition*)*nMarginals);
	muH=(double***) malloc(sizeof(double**)*nMarginals);
	alphaH=(double***) malloc(sizeof(double**)*nMarginals);
	betaH=(double***) malloc(sizeof(double**)*nMarginals);
	for(int i=0;i<nMarginals;i++) {
		betaH[i]=NULL;
	}
	costProvider=(THierarchicalCostFunctionProvider**) malloc(sizeof(THierarchicalCostFunctionProvider*)*nMarginals);
	weights=(double*) malloc(sizeof(double)*nMarginals);
	resH=(int**) malloc(sizeof(int*)*nMarginals);

}

void TMultiScaleSetupBarycenterContainer::cleanup() {
	if(HP!=NULL) {
		free(HP);
	}
	if(muH!=NULL) {
		free(muH);
	}
	if(alphaH!=NULL) {
		free(alphaH);
	}
	if(betaH!=NULL) {
		for(int i=0;i<nMarginals;i++) {
			if(betaH[i]!=NULL) {
				HPZ->signal_free_double(betaH[i], 0, HPZ->nLayers-1);
			}
		}
		free(betaH);
	}
	if(costProvider!=NULL) {
		free(costProvider);
	}
	if(weights!=NULL) {
		free(weights);
	}
	if(resH!=NULL) {
		free(resH);
	}
}

void TMultiScaleSetupBarycenterContainer::setMarginal(const int n, TMultiScaleSetupSingleBase &multiScaleSetup,
		const double weight) {
	HP[n]=multiScaleSetup.HP;
	muH[n]=multiScaleSetup.muH;
	alphaH[n]=multiScaleSetup.alphaH;
	weights[n]=weight;
	resH[n]=multiScaleSetup.resH;
}

void TMultiScaleSetupBarycenterContainer::setCostFunctionProvider(
		int n, THierarchicalCostFunctionProvider &costFunctionProvider) {
	costProvider[n]=&costFunctionProvider;
	costProvider[n]->alpha=alphaH[n];
	costProvider[n]->beta=betaH[n];
}

void TMultiScaleSetupBarycenterContainer::setCenterMarginal(TMultiScaleSetupSingleBase &multiScaleSetup) {
	HPZ=multiScaleSetup.HP;
	muZH=multiScaleSetup.muH;
	for(int i=0;i<nMarginals;i++) {
		betaH[i]=HPZ->signal_allocate_double(0,HPZ->nLayers-1);
	}
	resZH=multiScaleSetup.resH;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template class THierarchicalSearchList<THierarchicalNN::TCandidate>;

