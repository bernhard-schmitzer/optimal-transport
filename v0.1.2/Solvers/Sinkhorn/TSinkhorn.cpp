#include"TSinkhorn.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TKernelHandler_MultiCSR

template<class TMultiVarListHandlerType, class T>
TKernelHandler_MultiCSR_Base<TMultiVarListHandlerType,T>::TKernelHandler_MultiCSR_Base(int _dim, int *_res) {
	dim=_dim;
	varLists=(TMultiVarListHandlerType**) malloc(sizeof(TMultiVarListHandlerType*)*dim);

	int d;

	for(d=0;d<dim;d++) {
		varLists[d]=new TMultiVarListHandlerType(dim-1,_res[d]);
	}
}

template<class TMultiVarListHandlerType, class T>
TKernelHandler_MultiCSR_Base<TMultiVarListHandlerType,T>::~TKernelHandler_MultiCSR_Base() {
	int d;
	for(d=0;d<dim;d++) {
		delete varLists[d];
	}
	free(varLists);
}

template<class TMultiVarListHandlerType, class T>
void TKernelHandler_MultiCSR_Base<TMultiVarListHandlerType,T>::addVariables(int *pos, T value) {

	// add variable to all varLists
	int d,d2,offset;
	int *subPos;
	subPos=(int*) malloc(sizeof(int)*dim-1);
	for(d=0;d<dim;d++) {
		// add variable to varList[d]

		// copy elements of pos (except that at d) into subPos;
		offset=0;
		for(d2=0;d2<dim;d2++) {
			if(d2!=d) {
				subPos[offset]=pos[d2];
				offset++;
			}
		}

		varLists[d]->addToLine(pos[d], value, subPos);
	}
	free(subPos);
}

//template class TKernelHandler_MultiCSR_Base<TMultiVarListHandler<double>, double >;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TKernelHandler_CSR



TKernelHandler_CSR::TKernelHandler_CSR(int _dim, int *_res) {
	dim=_dim;
	varList=new TVarListHandler();

	varList->setupEmpty(_res[0]);
}

TKernelHandler_CSR::TKernelHandler_CSR(TVarListHandler *_varList) {
	varList=_varList;
	dim=2;
}
TKernelHandler_CSR::~TKernelHandler_CSR() {
	if(varList!=NULL) {
		delete varList;
	}
}

void TKernelHandler_CSR::addVariables(int *pos, __attribute__((unused))  double value) {
	varList->addToLine(pos[0],pos[1]);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// THierarchicalDualConstraintSearch

template<class TKernelHandlerType, class T>
THierarchicalDualConstraintSearchBase<TKernelHandlerType,T>::THierarchicalDualConstraintSearchBase(
		int _dim, THierarchicalPartition **_partitions, TMultiCostFunctionProvider *_costProvider)
		{
	dim=_dim;
	partitions=_partitions;
	costProvider=_costProvider;
}


template<class TKernelHandlerType, class T>
TKernelHandlerType* THierarchicalDualConstraintSearchBase<TKernelHandlerType,T>::SearchDualConstraints(
		T const slack, int layerBottom) {
	int d;
	int *pos;
	TKernelHandlerType *kernel;
	
	pos=(int*) malloc(sizeof(int)*dim);

	// create empty kernel object. use pos to store resolution of multi-marginal-array
	for(d=0;d<dim;d++) {
		pos[d]=partitions[d]->layers[layerBottom]->nCells;
	}
	kernel = new TKernelHandlerType(dim,pos);


	// reset pos to zero
	for(d=0;d<dim;d++) {
		pos[d]=0;
	}

	SearchDualConstraints_checkCell(slack, kernel, layerBottom, 0, pos);
	
	free(pos);
	
	return kernel;
	
}


template<class TKernelHandlerType, class T>
void THierarchicalDualConstraintSearchBase<TKernelHandlerType,T>::SearchDualConstraints_checkCell(
		T const slack, TKernelHandlerType * const kernel, int layerBottom,
		int layer, int *pos) {
	// on a given hierarchy level (specified by layer), on a given cell (specified by pos), check dual constraint
	// if found violated, launch refined search or add variable
	
	// start with slack equal to cost
	T slackValue=costProvider->getCostEff(layer, pos);


	// now check if slack is violated
	bool violation=(slackValue<=slack);

	// now, depending on violation, either refine search or add variables
	if(violation) {
		
		if(layer==layerBottom) {
			// add variable to kernel
			kernel->addVariables(pos, slackValue);
		} else {
			// start search on finer level
			SearchDualConstraints_checkChildren(slack, kernel,
				layerBottom, layer, pos);
		}
	}
}

template<class TKernelHandlerType, class T>
void THierarchicalDualConstraintSearchBase<TKernelHandlerType,T>::SearchDualConstraints_checkChildren(
		T const slack, TKernelHandlerType * const kernel,
		int layerBottom, int layerParent, int *posParent) {
	// recursively iterate over children of a given cell.
	
	int activeDim;
	int *posIndex;
	int *pos;

	// verbose

	pos=(int*) malloc(sizeof(int)*dim);
	posIndex=(int*) malloc(sizeof(int)*dim);

	// iterate recursively over all dimensions
	activeDim=0;
	posIndex[0]=0;
	while(activeDim>=0) {

		if (posIndex[activeDim]<(int) partitions[activeDim]->layers[layerParent]->nChildren[posParent[activeDim]]) {
			// if there are more children in current axis

			pos[activeDim]=partitions[activeDim]->layers[layerParent]->children[posParent[activeDim]]
			                                                                    [posIndex[activeDim]];
			if(activeDim<(int) dim-1) {
				// while not yet on last axis, initialize next axis
				activeDim+=1;
				posIndex[activeDim]=-1;
			} else {
				SearchDualConstraints_checkCell(slack, kernel,
									layerBottom, layerParent+1, pos);
			}
		} else {
			// done iterating at current axis, jump back one axis
			activeDim-=1;
		}
		if(activeDim>=0) {
			posIndex[activeDim]+=1;
		}
	}
	free(pos);
	free(posIndex);
}

template class THierarchicalDualConstraintSearchBase<TKernelHandler_Pos,double >;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TMultiVarListRefinementBase


template<class TKernelHandlerTypeIn, class TKernelHandlerTypeOut, class T>
TMultiVarListRefinementBase<TKernelHandlerTypeIn,TKernelHandlerTypeOut,T>::TMultiVarListRefinementBase(
		int _dim, THierarchicalPartition **_partitions, TMultiCostFunctionProvider *_costProvider)
		{
	dim=_dim;
	partitions=_partitions;
	costProvider=_costProvider;
}


template<class TKernelHandlerTypeIn, class TKernelHandlerTypeOut, class T>
TKernelHandlerTypeOut* TMultiVarListRefinementBase<TKernelHandlerTypeIn,TKernelHandlerTypeOut,T>::Refine(
		TKernelHandlerTypeIn *kernelTop, int layerTop) {
	/* kenrelTop is a kernel handler on the level layerTop.
		a kernel, one level finer, is produced. */
	int d;
	int *posParent;
	TKernelHandlerTypeOut *kernel;
	
	// posParent, pos are arrays to store the coordinate in product space
	posParent=(int*) malloc(sizeof(int)*dim);

	// create empty kernel object. use posParent to store resolution of multi-marginal-array
	for(d=0;d<dim;d++) {
		posParent[d]=partitions[d]->layers[layerTop+1]->nCells;
	}
	kernel = new TKernelHandlerTypeOut(dim,posParent);

	
	// iterate over entries of kernel
	typename TKernelHandlerTypeIn::iterator it=kernelTop->getIterator();

	while(kernelTop->iterate_pos(&it,posParent)) {
		// call refinement routine on that cell
		Refine_refineCell(kernel, layerTop, posParent);
	}
	
	free(posParent);
	
	return kernel;
	
}


template<class TKernelHandlerTypeIn, class TKernelHandlerTypeOut, class T>
void TMultiVarListRefinementBase<TKernelHandlerTypeIn,TKernelHandlerTypeOut,T>::Refine_refineCell(
		TKernelHandlerTypeOut * const kernel, int layerTop, int *posParent) {
	/* refine the cell specified by posParent on layerTop. for this iterate recursively over all
		lists of children (current list given by activeDim). On final dim actually add all children to var lists. */
	

	int activeDim;
	int *posIndex;
	int *pos;

	// verbose

	pos=(int*) malloc(sizeof(int)*dim);
	posIndex=(int*) malloc(sizeof(int)*dim);

	// iterate recursively over all dimensions
	activeDim=0;
	posIndex[0]=0;
	while(activeDim>=0) {

		if (posIndex[activeDim]<(int) partitions[activeDim]->layers[layerTop]->nChildren[posParent[activeDim]]) {
			// if there are more children in current axis

			pos[activeDim]=partitions[activeDim]->layers[layerTop]->children[posParent[activeDim]][posIndex[activeDim]];
			if(activeDim<(int) dim-1) {
				// while not yet on last axis, initialize next axis
				activeDim+=1;
				posIndex[activeDim]=-1;
			} else {
				kernel->addVariables(pos, costProvider->getCostEff(layerTop+1, pos));
			}
		} else {
			// done iterating at current axis, jump back one axis
			activeDim-=1;
		}
		if(activeDim>=0) {
			posIndex[activeDim]+=1;
		}
	}
	free(pos);
	free(posIndex);
}

template class TMultiVarListRefinementBase<TKernelHandler_CSR,TKernelHandler_Pos,double >;

template class TMultiVarListRefinementBase<TKernelHandler_Pos,TKernelHandler_Pos,double >;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TMultiVarListEvaluationBase


template<class TKernelHandlerTypeIn, class TKernelHandlerTypeOut, class T>
TMultiVarListEvaluationBase<TKernelHandlerTypeIn,TKernelHandlerTypeOut,T>::TMultiVarListEvaluationBase(
		int _dim, THierarchicalPartition **_partitions, TMultiCostFunctionProvider *_costProvider)
		{
	dim=_dim;
	partitions=_partitions;
	costProvider=_costProvider;
}


template<class TKernelHandlerTypeIn, class TKernelHandlerTypeOut, class T>
TKernelHandlerTypeOut* TMultiVarListEvaluationBase<TKernelHandlerTypeIn,TKernelHandlerTypeOut,T>::Evaluate(
		TKernelHandlerTypeIn *kernelOld, int layerTop) {
	/* varListTop is a MultiVarList on the level layerTop, with the first product space as reference space.
		a list of refined MultiVarLists, one level finer, is produced, one for each product space as reference space. */
	int d;
	int *pos;
	TKernelHandlerTypeOut *kernel;

	// posParent, pos are arrays to store the coordinate in product space
	pos=(int*) malloc(sizeof(int)*dim);

	// create empty kernel object. use pos to store resolution of multi-marginal-array
	for(d=0;d<dim;d++) {
		pos[d]=partitions[d]->layers[layerTop]->nCells;
	}
	kernel = new TKernelHandlerTypeOut(dim,pos);

	// iterate over entries of kernel
	typename TKernelHandlerTypeIn::iterator it=kernelOld->getIterator();

	while(kernelOld->iterate_pos(&it,pos)) {
		// add cell with new cost values to new kernel object
		kernel->addVariables(pos, costProvider->getCostEff(layerTop, pos));
	}

	free(pos);

	return kernel;

}



template class TMultiVarListEvaluationBase<TKernelHandler_CSR,TKernelHandler_Pos,double >;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TKLProjectionBase

template<class TKernelHandlerType, class T>
TKLProjectionBase<TKernelHandlerType,T>::TKLProjectionBase(TKernelHandlerType *_k, T **_scaling, int *_res) {
	k=_k;
	res=_res;
	scaling=_scaling;
	dim=k->dim;
}

template<class TKernelHandlerType, class T>
void TKLProjectionBase<TKernelHandlerType,T>::getConvolution(int axis, double *muOut) {
	int x,d;
	int *pos;
	T sig;
	
	pos=(int*) malloc(sizeof(int)*dim);

	// first set out array to zero
	for(x=0;x<res[axis];x++) {
		muOut[x]=0;
	}


	// iterate over all entries

	// iterate over entries of kernel
	typename TKernelHandlerType::iterator it=k->getIterator();

	while(k->iterate_both(&it,pos,&sig)) {
		for(d=0;d<dim;d++) {
			if(d!=axis) {
				// do not include scaling[axis] into calculations. safes some effort and
				// allows some things to be done in place
				sig*=scaling[d][pos[d]];
			}
		}
		muOut[pos[axis]]+=sig;
	}

	free(pos);
}


template<class TMultiVarListHandlerType, class T>
void TKLProjectionBase<TMultiVarListHandlerType,T>::project(int axis, double *mu) {
	int x;

	// store convolution to scaling field of axis
	getConvolution(axis,scaling[axis]);
	
	// transform convolution into new scaling
	for(x=0;x<res[axis];x++) {
		scaling[axis][x]=mu[x]/scaling[axis][x];
	}
}


template<class TMultiVarListHandlerType, class T>
void TKLProjectionBase<TMultiVarListHandlerType,T>::projectSoft_KL(int axis, double *mu, double q) {
	int x;

	// store convolution to scaling field of axis
	getConvolution(axis,scaling[axis]);

	// transform convolution into new scaling
	for(x=0;x<res[axis];x++) {
		scaling[axis][x]=pow(mu[x]/scaling[axis][x],q);
	}
}


template<class TMultiVarListHandlerType, class T>
void TKLProjectionBase<TMultiVarListHandlerType,T>::getMarginal(int axis, double *muOut) {
	int x;

	// store convolution to out
	getConvolution(axis,muOut);
	
	// transform out into marginal
	for(x=0;x<res[axis];x++) {
		muOut[x]=scaling[axis][x]*muOut[x];
	}
}

template class TKLProjectionBase<TKernelHandler_Pos,double >;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TDenseKernelGetterBase

template<class TKernelHandlerType, class T>
TDenseCostGetterBase<TKernelHandlerType,T>::TDenseCostGetterBase(
		int _dim, THierarchicalPartition **_partitions, TMultiCostFunctionProvider *_costProvider,
		int _layer)
		{
	dim=_dim;
	partitions=_partitions;
	costProvider=_costProvider;
	layer=_layer;
}


template<class TKernelHandlerType, class T>
TKernelHandlerType* TDenseCostGetterBase<TKernelHandlerType,T>::getCost() {

	int d;
	int *pos;
	TKernelHandlerType *kernel;

	pos=(int*) malloc(sizeof(int)*dim);

	// create empty kernel object. use pos to store resolution of multi-marginal-array
	for(d=0;d<dim;d++) {
		pos[d]=partitions[d]->layers[layer]->nCells;
	}
	kernel = new TKernelHandlerType(dim,pos);


	// iterate recursively over all dimensions
	int activeDim;
	activeDim=0;
	pos[0]=0;
	while(activeDim>=0) {

		if (pos[activeDim]<(int) partitions[activeDim]->layers[layer]->nCells) {
			// if there are more children in current axis

			if(activeDim<(int) dim-1) {
				// while not yet on last axis, initialize next axis
				activeDim+=1;
				pos[activeDim]=-1;
			} else {
				kernel->addVariables(pos, costProvider->getCostEff(layer, pos));
			}
		} else {
			// done iterating at current axis, jump back one axis
			activeDim-=1;
		}
		if(activeDim>=0) {
			pos[activeDim]+=1;
		}
	}
	free(pos);

	return kernel;
}

template class TDenseCostGetterBase<TKernelHandler_Pos, double >;
