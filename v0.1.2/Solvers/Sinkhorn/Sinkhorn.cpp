#include"Sinkhorn.h"


long int Setup_CostFunctionProvider_SquaredEuclidean(TInteger64Matrix *posAddr,
		TInteger64Matrix *radiiAddr, TInteger64Matrix *alphaAddr, int posDim, int lBottom, double weight) {

	int dim;
	double ***pos, ***radii, ***alpha;

	dim=posAddr->dimensions[0];
	if(alphaAddr->data[0]==0) {
		alpha=NULL;
	} else {
		alpha=(double***) alphaAddr->data;
	}
	pos=(double***) posAddr->data;
	radii=(double***) radiiAddr->data;

	TMultiCostFunctionProvider_SquaredEuclidean *costFunctionProvider;
	costFunctionProvider= new TMultiCostFunctionProvider_SquaredEuclidean(pos,radii,dim,posDim,lBottom,alpha,weight);

	return (long int) costFunctionProvider;
}

long int Setup_CostFunctionProvider_Color_SquaredEuclidean_RGB(TInteger64Matrix *posAddr,
		TInteger64Matrix *radiiAddr, TInteger64Matrix *alphaAddr, int posDim, int lBottom,
		int lTop, double colorWeight, int FR_mode, double FR_kappa) {

	int dim;
	double ***pos, ***radii, ***alpha;

	dim=posAddr->dimensions[0];
	if(alphaAddr->data[0]==0) {
		alpha=NULL;
	} else {
		alpha=(double***) alphaAddr->data;
	}
	pos=(double***) posAddr->data;
	radii=(double***) radiiAddr->data;

	TMultiCostFunctionProvider_Color_SquaredEuclidean_RGB *costFunctionProvider;
	costFunctionProvider= new TMultiCostFunctionProvider_Color_SquaredEuclidean_RGB(pos,radii,alpha,
			dim,posDim,lBottom,
			lTop,colorWeight,FR_mode, FR_kappa);

	return (long int) costFunctionProvider;
}

long int Setup_CostFunctionProvider_Color_SquaredEuclidean_HSV_HS(TInteger64Matrix *posAddr,
		TInteger64Matrix *radiiPosAddr, TInteger64Matrix *radiiValAddr, TInteger64Matrix *radiiHueAddr,
		TInteger64Matrix *alphaAddr, int posDim, int lBottom,
		int liftMode, double colorWeight, int FR_mode, double FR_kappa) {
	int dim;
	double ***pos, ***radiiPos, ***radiiHue, ***radiiVal, ***alpha;

	dim=posAddr->dimensions[0];
	if(alphaAddr->data[0]==0) {
		alpha=NULL;
	} else {
		alpha=(double***) alphaAddr->data;
	}
	pos=(double***) posAddr->data;
	radiiPos=(double***) radiiPosAddr->data;
	radiiHue=(double***) radiiHueAddr->data;
	radiiVal=(double***) radiiValAddr->data;

	TMultiCostFunctionProvider_Color_SquaredEuclidean_HSV_HS *costFunctionProvider;
	costFunctionProvider= new TMultiCostFunctionProvider_Color_SquaredEuclidean_HSV_HS(pos,
			radiiPos,radiiHue,radiiVal,alpha,
			dim,posDim,lBottom,
			liftMode,colorWeight,FR_mode, FR_kappa);

	return (long int) costFunctionProvider;
}

long int Setup_CostFunctionProvider_SquaredEuclideanWF(TInteger64Matrix *posAddr,
		TInteger64Matrix *radiiAddr, TInteger64Matrix *alphaAddr, int posDim, int lBottom, double delta, double cMax) {

	int dim;
	double ***pos, ***radii, ***alpha;

	dim=posAddr->dimensions[0];
	if(alphaAddr->data[0]==0) {
		alpha=NULL;
	} else {
		alpha=(double***) alphaAddr->data;
	}
	pos=(double***) posAddr->data;
	radii=(double***) radiiAddr->data;

	TMultiCostFunctionProvider_SquaredEuclideanWF *costFunctionProvider;
	costFunctionProvider= new TMultiCostFunctionProvider_SquaredEuclideanWF(pos,radii,dim,posDim,lBottom,alpha, delta,
			cMax);

	return (long int) costFunctionProvider;
}

long int Setup_CostFunctionProvider_SquaredEuclideanBarycenter(TInteger64Matrix *posAddr,
		TInteger64Matrix *radiiAddr, TInteger64Matrix *alphaAddr, TDoubleMatrix *lambda, int posDim, int lBottom) {

	int dim;
	double ***pos, ***radii, ***alpha;

	dim=posAddr->dimensions[0];
	alpha=(double***) alphaAddr->data;
	pos=(double***) posAddr->data;
	radii=(double***) radiiAddr->data;

	TMultiCostFunctionProvider_SquaredEuclideanBarycenter *costFunctionProvider;
	costFunctionProvider= new TMultiCostFunctionProvider_SquaredEuclideanBarycenter(pos,radii,
			lambda->data,dim,posDim,lBottom,alpha);

	return (long int) costFunctionProvider;
}


long int Setup_CostFunctionProvider_Coulomb(TInteger64Matrix *posAddr,
		TInteger64Matrix *radiiAddr, TInteger64Matrix *alphaAddr, int posDim, int lBottom, TDoubleMatrix *charges) {
	TMultiCostFunctionProvider_Coulomb *costFunctionProvider;

	int dim;
	double ***pos, ***radii, ***alpha;

	dim=posAddr->dimensions[0];
	alpha=(double***) alphaAddr->data;
	pos=(double***) posAddr->data;
	radii=(double***) radiiAddr->data;

	costFunctionProvider= new TMultiCostFunctionProvider_Coulomb(pos,radii,
			dim,posDim,lBottom,alpha,charges->data);


	return (long int) costFunctionProvider;
}


long int Setup_CostFunctionProvider_Reflector_Spherical(TInteger64Matrix *posAddr,
		TInteger64Matrix *radiiAddr, TInteger64Matrix *alphaAddr, int posDim, int lBottom) {


	int dim;
	double ***pos, ***radii, ***alpha;

	dim=posAddr->dimensions[0];
	if(alphaAddr->data[0]==0) {
		alpha=NULL;
	} else {
		alpha=(double***) alphaAddr->data;
	}
	pos=(double***) posAddr->data;
	radii=(double***) radiiAddr->data;

	TMultiCostFunctionProvider_Reflector_Spherical *costFunctionProvider;
	costFunctionProvider= new TMultiCostFunctionProvider_Reflector_Spherical(pos,radii,dim,posDim,lBottom,alpha);

	return (long int) costFunctionProvider;


}


long int Setup_CostFunctionProvider_Interpolator(long int coarseAddr, long int fineAddr,
		TInteger64Matrix *partitionAddr, double q, TInteger64Matrix *alphaAddr) {
	TMultiCostFunctionProvider_Interpolator *costFunctionProvider;
	TMultiCostFunctionProvider *coarse, *fine;
	THierarchicalPartition **partition;
	double ***alpha;

	coarse=(TMultiCostFunctionProvider*) coarseAddr;
	fine=(TMultiCostFunctionProvider*) fineAddr;
	partition=(THierarchicalPartition**) partitionAddr->data;
	alpha=(double***) alphaAddr->data;

	costFunctionProvider = new TMultiCostFunctionProvider_Interpolator(coarse,fine,partition,q,alpha, true);

	return (long int) costFunctionProvider;
}



int Check_DualConstraints_Pos(TInteger64Matrix *PartitionAddr,
		int lBottom,
		double slack,
		long int costFunctionProviderAddr,
		TInteger32Matrix *varListSpecs, TInteger64Matrix *varListPointer) {

	int dim;
	THierarchicalPartition **partitions;

	dim=PartitionAddr->dimensions[0];
	partitions=(THierarchicalPartition**) PartitionAddr->data;

	THierarchicalDualConstraintSearch_Pos *DualConstraintSearch;
	TMultiCostFunctionProvider *costFunctionProvider;

	costFunctionProvider=(TMultiCostFunctionProvider*) costFunctionProviderAddr;
	DualConstraintSearch = new THierarchicalDualConstraintSearch_Pos(dim,partitions,costFunctionProvider);

	TKernelHandler_Pos *result;
	result=DualConstraintSearch->SearchDualConstraints(slack,lBottom);

	varListSpecs->data[0]=result->total;
	varListSpecs->data[1]=result->dim;
	varListPointer->data[0]=(long int) result;

	delete DualConstraintSearch;
	return 0;
}


int Refine_VarList_CSR_Pos(TInteger64Matrix *PartitionAddr,
		int lTop,
		long int costFunctionProviderAddr,
		TInteger32Matrix *indices, TInteger32Matrix *indptr,
		TInteger32Matrix *varListSpecs, TInteger64Matrix *varListPointer) {

	int dim;
	THierarchicalPartition **partitions;

	dim=PartitionAddr->dimensions[0];
	partitions=(THierarchicalPartition**) PartitionAddr->data;

	TMultiCostFunctionProvider *costFunctionProvider;
	costFunctionProvider=(TMultiCostFunctionProvider*) costFunctionProviderAddr;

	TVarListHandler *varList;
	varList=new TVarListHandler();
	varList->fillFromCSRIndexList(indices->data,indptr->data,indptr->dimensions[0]-1,indices->dimensions[0]);
	TKernelHandler_CSR *kernel;
	kernel=new TKernelHandler_CSR(varList);

	TMultiVarListRefinement_CSR_Pos MultiVarListRefinement(dim,partitions,costFunctionProvider);

	TKernelHandler_Pos *kernelFine;
	kernelFine=MultiVarListRefinement.Refine(kernel,lTop);


	varListSpecs->data[0]=kernelFine->total;
	varListSpecs->data[1]=kernelFine->dim;
	varListPointer->data[0]=(long int) kernelFine;

	delete kernel;

	return 0;
}


int Refine_VarList_Pos_Pos(TInteger64Matrix *PartitionAddr,
		int lTop,
		long int costFunctionProviderAddr,
		long int kernelAddr,
		TInteger32Matrix *varListSpecs, TInteger64Matrix *varListPointer) {

	int dim;
	THierarchicalPartition **partitions;

	dim=PartitionAddr->dimensions[0];
	partitions=(THierarchicalPartition**) PartitionAddr->data;

	TMultiCostFunctionProvider *costFunctionProvider;
	costFunctionProvider=(TMultiCostFunctionProvider*) costFunctionProviderAddr;

	TKernelHandler_Pos *kernel;
	kernel=(TKernelHandler_Pos*) kernelAddr;

	TMultiVarListRefinement_Pos_Pos MultiVarListRefinement(dim,partitions,costFunctionProvider);

	TKernelHandler_Pos *kernelFine;
	kernelFine=MultiVarListRefinement.Refine(kernel,lTop);


	varListSpecs->data[0]=kernelFine->total;
	varListSpecs->data[1]=kernelFine->dim;
	varListPointer->data[0]=(long int) kernelFine;

	return 0;
}


int ReEvaluate_VarList_Pos(TInteger64Matrix *PartitionAddr,
		int lTop,
		long int costFunctionProviderAddr,
		TInteger32Matrix *indices, TInteger32Matrix *indptr,
		TInteger32Matrix *varListSpecs, TInteger64Matrix *varListPointer) {

	int dim;
	THierarchicalPartition **partitions;

	dim=PartitionAddr->dimensions[0];
	partitions=(THierarchicalPartition**) PartitionAddr->data;

	TMultiCostFunctionProvider *costFunctionProvider;
	costFunctionProvider=(TMultiCostFunctionProvider*) costFunctionProviderAddr;

	TVarListHandler *varList;
	varList=new TVarListHandler();
	varList->fillFromCSRIndexList(indices->data,indptr->data,indptr->dimensions[0]-1,indices->dimensions[0]);
	TKernelHandler_CSR *kernel;
	kernel=new TKernelHandler_CSR(varList);

	TMultiVarListEvaluation_CSR_Pos MultiVarListEvaluation(dim,partitions,costFunctionProvider);

	TKernelHandler_Pos *kernelNew;
	kernelNew=MultiVarListEvaluation.Evaluate(kernel,lTop);


	varListSpecs->data[0]=kernelNew->total;
	varListSpecs->data[1]=kernelNew->dim;
	varListPointer->data[0]=(long int) kernelNew;

	delete kernel;

	return 0;
}



int Iterate(long int kernelAddr, TInteger64Matrix *scalingAddr, TInteger64Matrix *muAddr,
		TInteger32Matrix *res, int n) {

	double **scaling, **mu;
	TKernelHandler_Pos *kernel;
	TKLProjection_Pos *projector;

	int dim;
	int d;


	kernel=(TKernelHandler_Pos*) kernelAddr;

	// some simply sanity checks
	if(scalingAddr->dimensions[0]!=muAddr->dimensions[0]) {
		return -1;
	}
	if(scalingAddr->dimensions[0]!=kernel->dim) {
		return -2;
	}
	if(scalingAddr->dimensions[0]!=res->dimensions[0]) {
		return -3;
	}

	// extract number of involved marginals
	dim=scalingAddr->dimensions[0];

	scaling=(double **) scalingAddr->data;
	mu=(double **) muAddr->data;

	projector=new TKLProjection_Pos(kernel,scaling,res->data);


	// do iterations
	int i;
	for(i=0;i<(int) n;i++) {
		for(d=0;d<dim;d++) {
			projector->project(d,mu[d]);
		}
	}

	// clean up
	delete projector;

	return 0;
}


int Tools_Collect_PosVarList(TDoubleMatrix *signal, TInteger32Matrix *indices,
		long int kernelPointer, int doDelete) {
	TKernelHandler_Pos *kernel;
	kernel=(TKernelHandler_Pos*) kernelPointer;

	// basic sanity check: do arrays have proper dimensions
	if(indices->depth!=2) {
		return-5;
	}
	if(kernel->total != indices->dimensions[0]) {
		return -2;
	}
	if(kernel->total != signal->dimensions[0]) {
		return -3;
	}
	if(kernel->dim != indices->dimensions[1]) {
		return -4;
	}

	kernel->writeToPosIndexList(signal->data, indices->data);

	if (doDelete>0) {
		delete kernel;
	}
	return 0;
}

int Tools_GetSpecs_PosVarList(long int kernelPointer, TInteger32Matrix *specs) {
	TKernelHandler_Pos *kernel;
	kernel=(TKernelHandler_Pos*) kernelPointer;

	specs->data[0]=kernel->total;
	specs->data[1]=kernel->dim;

	return 0;
}

int Tools_Export_PosVarList(TDoubleMatrix *signal, TInteger32Matrix *indices,
		TInteger64Matrix *pointer) {
	TKernelHandler_Pos *kernel;

	if(indices->depth!=2) {
		return -5;
	}
	if(signal->dimensions[0]!=indices->dimensions[0]) {
		return -4;
	}

	kernel=new TKernelHandler_Pos(indices->dimensions[1],NULL);
	kernel->fillFromPosIndexList(signal->data,indices->data,signal->dimensions[0]);
	pointer->data[0]=(long int) kernel;

	return 0;
}


int Tools_Delete_PosVarList(long int kernelPointer) {
	TKernelHandler_Pos *kernel;
	kernel=(TKernelHandler_Pos*) kernelPointer;
	delete kernel;
	return 0;
}


int Tools_Delete_CostFunctionProvider(long int CFPAddr) {
	TMultiCostFunctionProvider *CFP;
	CFP=(TMultiCostFunctionProvider*) CFPAddr;
	delete CFP;
	return 0;
}


// evaluate cost on a whole list of pos values. given by 2d array pos: one pos per row.
int Tools_Evaluate_CostFunctionProvider(long int CFPAddr, TInteger32Matrix *pos, TDoubleMatrix *c, int layer) {

	TMultiCostFunctionProvider *costFunctionProvider;

	costFunctionProvider=(TMultiCostFunctionProvider*) CFPAddr;

	int x;
	int dim;
	dim=pos->dimensions[1];
	// iterate over all rows
	for(x=0;x<pos->dimensions[0];x++) {
		c->data[x]=costFunctionProvider->getCost(layer,pos->data+x*dim);
	}
	return 0;

}



int Tools_GetMarginal(long int kernelAddr, TInteger64Matrix *scalingAddr, TDoubleMatrix *mu,
		TInteger32Matrix *res, int axis) {
	// axis specifies along which axis kernel is organized and for which axis the marginal is to be computed


	TKernelHandler_Pos *kernel;
	double **scaling;


	kernel=(TKernelHandler_Pos*) kernelAddr;

	scaling=(double **) scalingAddr->data;

	TKLProjection_Pos *projector;
	projector=new TKLProjection_Pos(kernel,scaling,res->data);
	projector->getMarginal(axis,mu->data);

	delete projector;

	return 0;
}



int Tools_GetDenseCosts_Pos(TInteger64Matrix *PartitionAddr, int layer,
		long int costFunctionProviderAddr,
		TInteger32Matrix *varListSpecs, TInteger64Matrix *varListPointer) {


	int dim;
	THierarchicalPartition **partitions;

	dim=PartitionAddr->dimensions[0];
	partitions=(THierarchicalPartition**) PartitionAddr->data;

	TMultiCostFunctionProvider *costFunctionProvider;
	costFunctionProvider=(TMultiCostFunctionProvider*) costFunctionProviderAddr;

	TDenseCostGetter_Pos *DenseCostGetter;
	DenseCostGetter = new TDenseCostGetter_Pos(dim,partitions, costFunctionProvider, layer);
	TKernelHandler_Pos *result;
	result=DenseCostGetter->getCost();

	varListSpecs->data[0]=result->total;
	varListSpecs->data[1]=result->dim;
	varListPointer->data[0]=(long int) result;

	delete DenseCostGetter;

	return 0;
}


int Kernel_Scale(long int kernelAddr, TInteger64Matrix *muAddr) {
	double **mu;
	TKernelHandler_Pos *kernel;

	kernel=(TKernelHandler_Pos*) kernelAddr;

	mu=(double **) muAddr->data;

	kernel->rescale(mu);

	return 0;

}

int Kernel_Exponentiate(long int kernelAddr, double eps) {
	TKernelHandler_Pos *kernel;

	kernel=(TKernelHandler_Pos*) kernelAddr;

	kernel->exponentiate(eps);

	return 0;

}
