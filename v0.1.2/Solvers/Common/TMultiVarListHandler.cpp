#include "TMultiVarListHandler.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TMultiVarListHandler


template <class T>
TMultiVarListHandler<T>::TMultiVarListHandler(int _dim) {
	res=0;
	total=0;
	dim=_dim;
	lenList=NULL;
	varList=NULL;
}

template <class T>
TMultiVarListHandler<T>::TMultiVarListHandler(int _dim, int _res) {
	res=0;
	total=0;
	dim=_dim;
	lenList=NULL;
	varList=NULL;
	setupEmpty(_res);
}


template <class T>
TMultiVarListHandler<T>::~TMultiVarListHandler() {
	clear();
}

template <class T>
void TMultiVarListHandler<T>::clear() {
	if(lenList!=NULL) {
			// go through varLists
			for(int i=0;i<res;i++) {
				// go through single var List
				for(int j=0;j<lenList->at(i);j++) {
					// free single coordinate list
					free(varList[i]->at(j));
				}
				delete varList[i];
				// free signal list
				delete signalList[i];
			}
			free(varList);
			delete lenList;
	}
	varList=NULL;
	lenList=NULL;
	res=0;
	total=0;
}

template <class T>
void TMultiVarListHandler<T>::setupEmpty(int _res) {
	res=_res;
	lenList=new vector<int>(res);
	varList=(vector<int*>**) malloc(sizeof(vector<int*>*)*res);
	signalList=(vector<T>**) malloc(sizeof(vector<T>*)*res);
	for(int x=0;x<res;x++) {
		varList[x]=new vector<int*>(0);
		signalList[x]=new vector<T>(0);
	}
}

template <class T>
void TMultiVarListHandler<T>::fillFromCSRIndexList(T *signal, int *indices, int *indptr, int _res, int _total) {
	setupEmpty(_res);
	total=_total;
	int rowLen,offset;
	int x,y,z;
	for(x=0;x<_res;x++) {
		
		rowLen=indptr[x+1]-indptr[x];
		offset=indptr[x];
		(*lenList)[x]=(int) rowLen;
		varList[x]->resize(rowLen);
		signalList[x]->resize(rowLen);
		for(y=0;y<rowLen;y++) {
			
			// allocate and assign coordinates
			(*(varList[x]))[y]=(int*) malloc(sizeof(int)*dim);


			for(z=0;z<(int) dim;z++) {
				((*(varList[x]))[y])[z]=indices[(offset+y)*dim+z];
				
			}
			// assign signal
			(*(signalList[x]))[y]=signal[offset+y];
			
		}
	}
}

template <class T>
void TMultiVarListHandler<T>::writeToCSRIndexList(T *signal, int *indices, int *indptr) {
	// write varList data to CSR (as in scipy.sparse.csr_matrix) list. indices must have size total, indptr must have size res+1
	int x,yIndex,z,offset;
	indptr[0]=0;
	offset=0;
	for(x=0;x<res;x++) {
		for(yIndex=0;yIndex<(*lenList)[x];yIndex++) {
			for(z=0;z<dim;z++) {
				// assign coordinates
				indices[offset*dim+z]=((*(varList[x]))[yIndex])[z];
				// assign signal
				signal[offset]=(*(signalList[x]))[yIndex];
			}
			offset++;
		}
		indptr[x+1]=offset;
	}
	
}
template <class T>
void TMultiVarListHandler<T>::addToLine(int x, T signal, int *yCandidate) {
	int yIndex,z;
	// assign coordinates
	varList[x]->push_back((int*) malloc(sizeof(int)*dim));
	yIndex=lenList->at(x);
	for(z=0;z<dim;z++) {
		((*(varList[x]))[yIndex])[z]=yCandidate[z];
	}
	// assign signal
	signalList[x]->push_back(signal);
	lenList->at(x)++;
	total++;
}

template class TMultiVarListHandler<int>;
template class TMultiVarListHandler<double>;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TMultiCostFunctionProvider


TMultiCostFunctionProvider::TMultiCostFunctionProvider(double ***_pos, double ***_radii, int _dim, int _posDim, int _layerBottom,
		double *** _alpha) {
	pos=_pos;
	radii=_radii;
	dim=_dim;
	posDim=_posDim;
	layerBottom=_layerBottom;
	alpha=_alpha;
}

TMultiCostFunctionProvider::~TMultiCostFunctionProvider() {
}



///////////////////////////////////////
// Euclidean 2 marginals


TMultiCostFunctionProvider_SquaredEuclidean::TMultiCostFunctionProvider_SquaredEuclidean(
		double ***_pos, double ***_radii, int _dim, int _posDim, int _layerBottom, double ***_alpha, double _weight) :
		TMultiCostFunctionProvider(_pos, _radii, _dim, _posDim, _layerBottom, _alpha) {
		
	weight=_weight;
}

double TMultiCostFunctionProvider_SquaredEuclidean::getCost(int layer, int *x) {
	double result;
	
	// compute squared Euclidean distance between representatives, fixed to dim=2
	result=EUCL_lincombSqr(pos[0][layer]+(x[0]*posDim), pos[1][layer]+(x[1]*posDim), 1, -1, posDim);

	if(layer<layerBottom) {
		// if not at finest layer, need to compute lower bound
		result=sqrt(result)-radii[0][layer][x[0]]-radii[1][layer][x[1]];
		if(result>=0) {
			result=pow(result,2);
		} else {
			result=0;
		}
	}

	return weight*result;

}


///////////////////////////////////////
// Euclidean + Fisher Rao, 2 marginals


TMultiCostFunctionProvider_SquaredEuclideanWF::TMultiCostFunctionProvider_SquaredEuclideanWF(
		double ***_pos, double ***_radii, int _dim, int _posDim, int _layerBottom, double ***_alpha, double _delta, double _cMax) :
		TMultiCostFunctionProvider(_pos, _radii, _dim, _posDim, _layerBottom, _alpha) {
		
		delta=_delta;
		cMax=_cMax;
		prefactor=2*delta*delta;
}

double TMultiCostFunctionProvider_SquaredEuclideanWF::getCost(int layer, int *x) {
	double result;
	
	// compute Euclidean distance between representatives, fixed to dim=2
	result=EUCL_lincombSqr(pos[0][layer]+(x[0]*posDim), pos[1][layer]+(x[1]*posDim), 1, -1, posDim);
	result=sqrt(result);

	if(layer<layerBottom) {
		// if not at finest layer, need to compute lower bound
		result=result-radii[0][layer][x[0]]-radii[1][layer][x[1]];
		if(result<0) {
			result=0;
		}
	}
	
	// now post-processing for OT-FR cost
	
	if(result>M_PI*delta) {
		return prefactor*cMax;
	} else {
		return prefactor*min(cMax,-2*log(cos(result/(2*delta))));
	}

}

///////////////////////////////////
// Barycenter

TMultiCostFunctionProvider_SquaredEuclideanBarycenter::TMultiCostFunctionProvider_SquaredEuclideanBarycenter(
		double ***_pos, double ***_radii, double *_lambda, int _dim, int _posDim, int _layerBottom, double ***_alpha) :
		TMultiCostFunctionProvider(_pos, _radii, _dim, _posDim, _layerBottom, _alpha) {
	lambda=_lambda;
}

double TMultiCostFunctionProvider_SquaredEuclideanBarycenter::getCost(int layer, int *x) {
	int i,j;
	double result;
	double radSum;
	double xAbs;

	// compute exact Barycenter cost (only the inner product part, only between different indices, only one triangle of matrix)
	result=0;
	for(i=0;i<dim;i++) {
		for(j=0;j<i;j++) {
			result-=lambda[i]*lambda[j]*EUCL_innerProduct(pos[i][layer]+(x[i]*posDim), pos[j][layer]+(x[j]*posDim),posDim);
		}
	}
	
	if(layer<layerBottom) {
		// if not at finest layer, need to compute lower bound
	
		// r_i r_j contribution
		for(i=0;i<dim;i++) {
			for(j=0;j<i;j++) {
				result-=lambda[i]*lambda[j]*radii[j][layer][x[j]]*radii[j][layer][x[j]];
			}
		}
	
		// |x_i| r_j contribution
		for(i=0;i<dim;i++) {
			radSum=0;
			for(j=0;j<dim;j++) {
				if(j!=i) {
					radSum+=lambda[j]*radii[j][layer][x[j]];
				}
			}
			xAbs=sqrt(EUCL_innerProduct(pos[i][layer]+(x[i]*posDim),pos[i][layer]+(x[i]*posDim), posDim));
			result-=lambda[i]*xAbs*radSum;
		}
	
	}

	return result;
}


///////////////////////////////////
// Coulomb

TMultiCostFunctionProvider_Coulomb::TMultiCostFunctionProvider_Coulomb(double ***_pos, double ***_radii,
		int _dim, int _posDim, int _layerBottom, double *** _alpha,
		double *_charges) :
		TMultiCostFunctionProvider(_pos, _radii, _dim, _posDim, _layerBottom, _alpha) {
	charges=_charges;
}

double TMultiCostFunctionProvider_Coulomb::getCost(int layer, int *x) {
	int i,j;
	double result;
	double distance;
	
	result=0;
	
	for(i=0;i<dim;i++) {
		//cout << i << endl;
		for(j=i+1;j<dim;j++) {
			//cout << "\t" << j << endl;
			// compute exact distance between charges
			distance=sqrt(EUCL_lincombSqr(pos[i][layer]+(x[i]*posDim), pos[j][layer]+(x[j]*posDim), 1, -1, posDim));
			//cout << "\t\t" << distance << endl;

			if(layer<layerBottom) {
				// depending on sign of charges, get upper or lower bound of distance
				if(charges[i]*charges[j]>0) {
					//cout << "\t\t+" << endl;
					// if repulsive, make radius maximal
					distance+=radii[i][layer][x[i]];
					distance+=radii[j][layer][x[j]];
				} else {
					//cout << "\t\t-" << endl;
					distance-=radii[i][layer][x[i]];
					distance-=radii[j][layer][x[j]];
				}
				//cout << "\t\t" << distance << endl;
			}
			
			if(distance>0) {
				result+=charges[i]*charges[j]/distance;
			} else {
				result+=charges[i]*charges[j]*COST_INFINITY;
			}
				
			
		}
	}
	
	return result;
}

///////////////////////////////////
// Spherical Reflector


TMultiCostFunctionProvider_Reflector_Spherical::TMultiCostFunctionProvider_Reflector_Spherical(
		double ***_pos, double ***_radii, int _dim, int _posDim, int _layerBottom, double ***_alpha) :
		TMultiCostFunctionProvider(_pos, _radii, _dim, _posDim, _layerBottom, _alpha) {
}

double TMultiCostFunctionProvider_Reflector_Spherical::getCost(int layer, int *x) {
	double result;
	
	// compute cos phi
	result=EUCL_innerProduct(pos[0][layer]+(x[0]*posDim), pos[1][layer]+(x[1]*posDim), posDim);
	if(result>1.) {
		result=1.;
	}
	if(result<-1.) {
		result=-1.;
	}

	if(layer<layerBottom) {
		// if not at finest layer, need to compute upper bound on phi

		// get phi
		result=acos(result);
		// add radii to get upper bound
		result=result+radii[0][layer][x[0]]+radii[1][layer][x[1]];
		// truncate if phi bound is larger than pi
		if(result>M_PI) {
			result=M_PI;
		}
		// recompuse cos phi
		result=cos(result);		
	}

	// evaluate actual cost function
	result=-log(1.-result);

	return result;

}


///////////////////////////////////
// Interpolator

TMultiCostFunctionProvider_Interpolator::TMultiCostFunctionProvider_Interpolator(
		TMultiCostFunctionProvider *_coarse, TMultiCostFunctionProvider *_fine,
		THierarchicalPartition **_partition, double _q, double *** _alpha, bool _destroyChildren) :
		TMultiCostFunctionProvider(NULL, NULL, _fine->dim, _fine->posDim, _fine->layerBottom, _alpha) {
	coarse=_coarse;
	fine=_fine;
	partition=_partition;
	q=_q;
	destroyChildren=_destroyChildren;
}

TMultiCostFunctionProvider_Interpolator::~TMultiCostFunctionProvider_Interpolator() {
	if (destroyChildren) {
		delete coarse;
		delete fine;
	}
}


double TMultiCostFunctionProvider_Interpolator::getCost(int layer, int *x) {
	// if on top layer, so no coarser rep can be given, return basically -inf
	if (layer==0) {
		return -DBL_INFINITY;
	}
	
	// determine pos of parent element
	int d;
	int *posParent;
	posParent=(int*) malloc(sizeof(int)*dim);
	for(d=0;d<dim;d++) {
		posParent[d]=partition[d]->layers[layer]->parent[x[d]];
	}
	
	double result;
	result=q*fine->getCost(layer,x)+(1-q)*coarse->getCost(layer-1,posParent);
	free(posParent);
	return result;
}


/////////////////////////////////////////////////////////////////
// Color Transport

///////////////////////////////////////
/* Squared Euclidean + RGB cyclic */


TMultiCostFunctionProvider_Color_SquaredEuclidean_RGB::TMultiCostFunctionProvider_Color_SquaredEuclidean_RGB(
		double ***_pos, double ***_radii, double *** _alpha,
		int _dim, int _posDim, int _layerBottom,
		int _lTop, double _colorWeight, bool _FR_mode, double _FR_delta) :
		TMultiCostFunctionProvider(_pos, _radii, _dim, _posDim, _layerBottom, _alpha) {
		
	colorWeight=_colorWeight;
	lTop=_lTop;
	FR_mode=_FR_mode;
	FR_delta=_FR_delta;
	FR_prefactor=2*FR_delta*FR_delta;

}

double TMultiCostFunctionProvider_Color_SquaredEuclidean_RGB::getCost(int layer, int *x) {
	double distEucl,distColor;
	double result;
	
	if(layer<lTop) {
		// if too fine for partition to have properly separated colors, simply return 0.
		return 0;
	}
	
	// compute Euclidean distance between spatial coordinates
	distEucl=EUCL_lincombSqr(pos[0][layer]+(x[0]*posDim), pos[1][layer]+(x[1]*posDim), 1, -1, posDim-1);

	if(layer<layerBottom) {
		// if not at finest layer, need to compute lower bound
		distEucl=sqrt(distEucl)-radii[0][layer][x[0]]-radii[1][layer][x[1]];
		if(distEucl>=0) {
			distEucl=pow(distEucl,2);
		} else {
			distEucl=0;
		}
	}
	
	// check distance in color dimension
	if(abs(pos[0][layer][x[0]*posDim+(posDim-1)]-pos[1][layer][x[1]*posDim+(posDim-1)])<1E-10) {
		// if sufficiently close, assume is same channel
		distColor=0;
	} else {
		distColor=colorWeight*colorWeight;
	}

	result=distEucl+distColor;
	if(FR_mode) {
		return getWFLogCost(sqrt(result), FR_delta, FR_prefactor);
	} else {
		return result;
	}

}


///////////////////////////////////////
/* Squared Euclidean + HSV_HS */


TMultiCostFunctionProvider_Color_SquaredEuclidean_HSV_HS::TMultiCostFunctionProvider_Color_SquaredEuclidean_HSV_HS(
		double ***_pos,
		double ***_radiiSpace, double ***_radiiHue, double ***_radiiVal,
		double *** _alpha,
		int _dim, int _posDim, int _layerBottom,
		int _liftMode,
		double _colorWeight, bool _FR_mode, double _FR_delta) :
		TMultiCostFunctionProvider(_pos, _radiiSpace, _dim, _posDim, _layerBottom, _alpha) {
		
	liftMode=_liftMode;
	colorWeight=_colorWeight;
	FR_mode=_FR_mode;
	FR_delta=_FR_delta;
	FR_prefactor=2*FR_delta*FR_delta;
	
	radiiHue=_radiiHue;
	radiiVal=_radiiVal;

}

double TMultiCostFunctionProvider_Color_SquaredEuclidean_HSV_HS::getCost(int layer, int *x) {
	double distEucl,distVal,distHue;
	double result;
	
	// compute Euclidean distance between spatial coordinates
	distEucl=EUCL_lincombSqr(pos[0][layer]+(x[0]*posDim), pos[1][layer]+(x[1]*posDim), 1, -1, posDim-2);

	if(layer<layerBottom) {
		// if not at finest layer, need to compute lower bound
		distEucl=sqrt(distEucl)-radii[0][layer][x[0]]-radii[1][layer][x[1]];
		if(distEucl>=0) {
			distEucl=pow(distEucl,2);
		} else {
			distEucl=0;
		}
	}
	
	// distance in val dimension
	distVal=EUCL_lincombSqr(pos[0][layer]+(x[0]*posDim)+(posDim-1), pos[1][layer]+(x[1]*posDim)+(posDim-1), 1, -1, 1);
	if(layer<layerBottom) {
		// if not at finest layer, need to compute lower bound
		distVal=sqrt(distVal)-radiiVal[0][layer][x[0]]-radiiVal[1][layer][x[1]];
		if(distVal>=0) {
			distVal=pow(distVal,2);
		} else {
			distVal=0;
		}
	}
	

	// S1-type distance in hue dimension
	distHue=EUCL_lincombSqr(pos[0][layer]+(x[0]*posDim)+(posDim-2), pos[1][layer]+(x[1]*posDim)+(posDim-2), 1, -1, 1);
	// S1-minimization
	distHue=sqrt(distHue);
	distHue=min(distHue,1-distHue);
	if(layer<layerBottom) {
		// if not at finest layer, need to compute lower bound
		distHue=distHue-radiiHue[0][layer][x[0]]-radiiHue[1][layer][x[1]];
		if(distHue<0) {
			distHue=0;
		}
	}
	// re-square
	distHue=pow(distHue,2);

	

	result=distEucl+colorWeight*colorWeight*(distVal+distHue);
	if(FR_mode) {
		return getWFLogCost(sqrt(result), FR_delta, FR_prefactor);
	} else {
		return result;
	}

}

