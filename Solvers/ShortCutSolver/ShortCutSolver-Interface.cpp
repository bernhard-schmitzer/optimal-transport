#include"ShortCutSolver-Interface.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Cost Function Provider Methods */


int Setup_CostFunctionProvider_SqrEuclidean(
		TDoubleMatrix *xPos, TDoubleMatrix *yPos,
		TInteger64Matrix *Pointer) {

	int xres,yres;
	int dim;

	xres=xPos->dimensions[0];
	yres=yPos->dimensions[0];
	dim=xPos->dimensions[1];


	TCostFunctionProvider_Dynamic *costFunctionProvider;
	costFunctionProvider=new TCostFunctionProvider_Dynamic(xres, yres, xPos->data, yPos->data, dim);

	// return pointers
	Pointer->data[0]=(long int) costFunctionProvider;

	return 0;
}


int Setup_CostFunctionProvider_Torus(
		TDoubleMatrix *xPos, TDoubleMatrix *yPos, TDoubleMatrix *radius, int torusDim,
		TInteger64Matrix *Pointer) {

	int xres,yres;
	int dim;

	xres=xPos->dimensions[0];
	yres=yPos->dimensions[0];
	dim=xPos->dimensions[1];


	TCostFunctionProvider_Dynamic_Torus *costFunctionProvider;
	costFunctionProvider=new TCostFunctionProvider_Dynamic_Torus(xres, yres, xPos->data, yPos->data, dim,
			radius->data, torusDim);

	// return pointers
	Pointer->data[0]=(long int) costFunctionProvider;

	return 0;
}


/* Interaction with CostFunctionProvider */
 int CostFunctionProvider_Evaluate(long int costFunctionProviderAddr,
		 	TInteger32Matrix *VLXindices, TInteger32Matrix *VLXindptr, TDoubleMatrix *c) {

	TVarListHandler *xVars;
	xVars=SetupXVars(VLXindices->data, VLXindptr->data, VLXindptr->dimensions[0]-1, VLXindices->dimensions[0]);


	TCostFunctionProviderBase *cProvider;
	cProvider=(TCostFunctionProviderBase*) costFunctionProviderAddr;

	double *cReturn;
	cReturn=cProvider->getC(xVars);

	for(int i=0;i<VLXindices->dimensions[0];i++) {
		c->data[i]=cReturn[i];
	}

	delete xVars;

	if(cProvider->free_c()) {
		free(cReturn);
	}
	return 0;
 }


 int CostFunctionProvider_GetDenseCost(long int costFunctionProviderAddr,
 		 	TDoubleMatrix *c) {

	TCostFunctionProvider_Dynamic *cProvider;
	cProvider=(TCostFunctionProvider_Dynamic*) costFunctionProviderAddr;

	cProvider->getCDense(c->data);
	return 0;
}


 int CostFunctionProvider_GetRes(long int costFunctionProviderAddr,
 		 	TInteger32Matrix *res) {

	 TCostFunctionProvider_Dynamic *cProvider;
	 cProvider=(TCostFunctionProvider_Dynamic*) costFunctionProviderAddr;
	 res->data[0]=cProvider->xres;
	 res->data[1]=cProvider->yres;
	 return 0;
 }


 int CostFunctionProvider_Delete(long int costFunctionProviderAddr) {

 	 TCostFunctionProvider_Dynamic *cProvider;
 	 cProvider=(TCostFunctionProvider_Dynamic*) costFunctionProviderAddr;
 	 delete cProvider;
 	 return 0;
  }


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Coupling Handler Methods */


int Setup_CouplingHandler_SemiDense(TDoubleMatrix* c, TDoubleMatrix* mu,
		TInteger32Matrix* VLXindices, TInteger32Matrix* VLXindptr,
		TInteger64Matrix* Pointer) {

	int xres,yres;

	xres=c->dimensions[0];
	yres=c->dimensions[1];

	TVarListHandler *xVars;
	xVars=SetupXVars(VLXindices->data, VLXindptr->data, xres, VLXindices->dimensions[0]);

	// coupling handler setup
	TCouplingHandlerSemiDense *couplingHandler;
	couplingHandler = new TCouplingHandlerSemiDense(xres, yres, c->data, mu->data, xVars);
	TCouplingHandlerExt<TCouplingHandlerSemiDense> *couplingHandlerInterface;
	couplingHandlerInterface = new TCouplingHandlerExt<TCouplingHandlerSemiDense>(couplingHandler,true);

	// return pointers
	Pointer->data[0]=(long int) couplingHandlerInterface;

	return 0;
}

int Setup_CouplingHandler_Sparse_fullC(TDoubleMatrix *c,
		TInteger32Matrix *VLXindices, TInteger32Matrix *VLXindptr,
		TInteger64Matrix *Pointer) {
	int xres,yres;

	xres=c->dimensions[0];
	yres=c->dimensions[1];

	TVarListHandler *xVars;
	xVars=SetupXVars(VLXindices->data, VLXindptr->data, xres, VLXindices->dimensions[0]);

	// cost function provider
	TCostFunctionProviderFullArray *costFunctionProvider;
	costFunctionProvider=new TCostFunctionProviderFullArray(xres,yres,c->data);

	// coupling handler setup
	TCouplingHandlerSparse *couplingHandler;
	couplingHandler = new TCouplingHandlerSparse(xres,yres, costFunctionProvider, xVars);
	TCouplingHandlerExt<TCouplingHandlerSparse> *couplingHandlerInterface;
	couplingHandlerInterface = new TCouplingHandlerExt<TCouplingHandlerSparse>(couplingHandler,true);


	// return pointers
	Pointer->data[0]=(long int) couplingHandlerInterface;

	return 0;
}



int Setup_CouplingHandler_Sparse_dynamicC(
		TInteger32Matrix *VLXindices, TInteger32Matrix *VLXindptr,
		int xres, int yres, long int costFunctionProviderAddr,
		TInteger64Matrix *Pointer) {

	TVarListHandler *xVars;
	xVars=SetupXVars(VLXindices->data, VLXindptr->data, xres, VLXindices->dimensions[0]);

	TCostFunctionProviderBase *costFunctionProvider;
	costFunctionProvider=(TCostFunctionProviderBase*) costFunctionProviderAddr;

	// coupling handler setup
	TCouplingHandlerSparse *couplingHandler;
	couplingHandler = new TCouplingHandlerSparse(xres,yres, costFunctionProvider, xVars);
	TCouplingHandlerExt<TCouplingHandlerSparse> *couplingHandlerInterface;
	couplingHandlerInterface = new TCouplingHandlerExt<TCouplingHandlerSparse>(couplingHandler,true);

	// return pointers
	Pointer->data[0]=(long int) couplingHandlerInterface;

	return 0;
}

/* Interaction with sparse coupling handler */

int CouplingHandler_Sparse_GetMu_Request(
		long int couplingHandlerAddr, TInteger32Matrix *MuSpecs) {
	TCouplingHandlerExt<TCouplingHandlerSparse> *couplingHandlerInterface;
	couplingHandlerInterface =(TCouplingHandlerExt<TCouplingHandlerSparse>*) couplingHandlerAddr;
	MuSpecs->data[0]=couplingHandlerInterface->couplingHandler->xres;
	MuSpecs->data[1]=couplingHandlerInterface->couplingHandler->total;
	return 0;
}

int CouplingHandler_Sparse_GetMu_Collect(
		long int couplingHandlerAddr, TDoubleMatrix *data, TInteger32Matrix *indices, TInteger32Matrix *indptr) {
	TCouplingHandlerExt<TCouplingHandlerSparse> *couplingHandlerInterface;
	couplingHandlerInterface =(TCouplingHandlerExt<TCouplingHandlerSparse>*) couplingHandlerAddr;

	couplingHandlerInterface->couplingHandler->xVars->writeToCSRIndexList(indices->data,indptr->data);

	int x;
	for(x=0;x<couplingHandlerInterface->couplingHandler->total;x++) {
		data->data[x]=couplingHandlerInterface->couplingHandler->mu[x];
	}

	return 0;
}


int CouplingHandler_Sparse_GetSupport_Request(long int couplingHandlerAddr, TInteger64Matrix *SupportPointer,
		TInteger32Matrix *SupportSpecs) {
	TCouplingHandlerExt<TCouplingHandlerSparse> *couplingHandlerInterface;
	couplingHandlerInterface =(TCouplingHandlerExt<TCouplingHandlerSparse>*) couplingHandlerAddr;

	TVarListSignal<double> *support;
	support=couplingHandlerInterface->getSupportSignal();

	SupportPointer->data[0]=(long int) support;
	SupportSpecs->data[0]=support->varList->res;
	SupportSpecs->data[1]=support->varList->total;

	return 0;
}

int CouplingHandler_Sparse_GetSupport_Collect(
		long int supportAddr,
		TDoubleMatrix *data, TInteger32Matrix *indices, TInteger32Matrix *indptr) {

	TVarListSignal<double> *support;
	support=(TVarListSignal<double>*) supportAddr;
	support->varList->writeToCSRIndexList(indices->data,indptr->data);

	int x;
	for(x=0;x<support->varList->total;x++) {
		data->data[x]=support->signal[x];
	}

	free(support->signal);
	delete support->varList;
	delete support;

	return 0;
}


int CouplingHandler_Sparse_GetCost(
		long int couplingHandlerAddr,
		TDoubleMatrix *data,
		int total) {

	TCouplingHandlerExt<TCouplingHandlerSparse> *couplingHandlerInterface;
	couplingHandlerInterface =(TCouplingHandlerExt<TCouplingHandlerSparse>*) couplingHandlerAddr;
	int yIndex;
	for(yIndex=0;yIndex<total;yIndex++) {
		data->data[yIndex]=couplingHandlerInterface->couplingHandler->c[yIndex];
	}
	return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Shielding Methods */

int Setup_Shielding_Grid(TInteger32Matrix *xDims, TInteger32Matrix *yDims, TInteger64Matrix *Pointer) {
	TShieldGeneratorGrid_SqrEuclidean *shieldGenerator;
	shieldGenerator= new TShieldGeneratorGrid_SqrEuclidean(xDims->dimensions[0],xDims->data,yDims->data);

	// return pointers
	Pointer->data[0]=(long int) shieldGenerator;
	return 0;
}


int Setup_Shielding_Padding(TInteger32Matrix *xDims, TInteger32Matrix *yDims, TInteger64Matrix *Pointer, int width) {
	TShieldGeneratorGrid_Padding *shieldGenerator;
	shieldGenerator= new TShieldGeneratorGrid_Padding(xDims->dimensions[0],xDims->data,yDims->data, width);

	// return pointers
	Pointer->data[0]=(long int) shieldGenerator;
	return 0;
}


int Setup_Shielding_Tree(
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		int lBottom, int lTop,
		TDoubleMatrix *xpos,
		long int yPartitionAddr, TInteger64Matrix *yPosAddr, TInteger64Matrix *yRadiiAddr,
		TInteger64Matrix *Pointer) {

	//int xres,yres,dim,msg;
	int xres,dim;

	xres=xpos->dimensions[0];
	dim=xpos->dimensions[1];

	// shielding generator setup
	TVarListHandler *xNeighbours;
	xNeighbours=new TVarListHandler();
	xNeighbours->fillFromCSRIndexList(NeighXindices->data, NeighXindptr->data, xres, NeighXindices->dimensions[0]);


	THierarchicalPartition *yPartition=(THierarchicalPartition*) yPartitionAddr;
	double **yPos= (double**) yPosAddr->data;
	double **yRadii = (double**) yRadiiAddr->data;

	TShieldGeneratorTree_SqrEuclidean *shieldGenerator;
	shieldGenerator= new TShieldGeneratorTree_SqrEuclidean(dim, yPartition, yPos, yRadii,
			lBottom, lTop,
			xpos->data, xNeighbours);


	Pointer->data[0]=(long int) shieldGenerator;

	return 0;
}

int Setup_Shielding_Tree_Torus(
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		int lBottom, int lTop,
		TDoubleMatrix *xpos,
		long int yPartitionAddr, TInteger64Matrix *yPosAddr, TInteger64Matrix *yRadiiAddr,
		TInteger64Matrix *yTorusRadiiAddr,
		TDoubleMatrix *torusRadii, int torusDim,
		TInteger64Matrix *Pointer) {

	//int xres,yres,dim,msg;
	int xres,dim;

	xres=xpos->dimensions[0];
	dim=xpos->dimensions[1];

	// shielding generator setup
	TVarListHandler *xNeighbours;
	xNeighbours=new TVarListHandler();
	xNeighbours->fillFromCSRIndexList(NeighXindices->data, NeighXindptr->data, xres, NeighXindices->dimensions[0]);


	THierarchicalPartition *yPartition=(THierarchicalPartition*) yPartitionAddr;
	double **yPos= (double**) yPosAddr->data;
	double **yRadii = (double**) yRadiiAddr->data;
	double ***yTorusRadii = (double***) yTorusRadiiAddr->data;

	TShieldGeneratorTree_Torus *shieldGenerator;
	shieldGenerator= new TShieldGeneratorTree_Torus(dim, yPartition, yPos, yRadii, yTorusRadii,
			lBottom, lTop,
			xpos->data, xNeighbours,
			torusRadii->data, torusDim);


	Pointer->data[0]=(long int) shieldGenerator;

	return 0;
}


int Setup_Shielding_TreeNoise(
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		int lBottom, int lTop,
		TDoubleMatrix *xpos,
		long int yPartitionAddr, TInteger64Matrix *yPosAddr, TInteger64Matrix *yRadiiAddr,
		TInteger64Matrix *cAddr, double eta, double lambda,
		TInteger64Matrix *Pointer) {

	//int xres,yres,dim,msg;
	int xres,dim;

	xres=xpos->dimensions[0];
	dim=xpos->dimensions[1];

	// shielding generator setup
	TVarListHandler *xNeighbours;
	xNeighbours=new TVarListHandler();
	xNeighbours->fillFromCSRIndexList(NeighXindices->data, NeighXindptr->data, xres, NeighXindices->dimensions[0]);


	THierarchicalPartition *yPartition=(THierarchicalPartition*) yPartitionAddr;
	double **yPos= (double**) yPosAddr->data;
	double **yRadii = (double**) yRadiiAddr->data;
	double **c = (double**) cAddr->data;

	TShieldGeneratorTree_SqrEuclideanNoise *shieldGenerator;
	shieldGenerator= new TShieldGeneratorTree_SqrEuclideanNoise(dim, yPartition, yPos, yRadii,
			lBottom, lTop,
			xpos->data, xNeighbours,
			c, eta, lambda);


	Pointer->data[0]=(long int) shieldGenerator;

	return 0;
}



int Setup_Shielding_TreePEucl(
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		int lBottom, int lTop,
		TDoubleMatrix *xpos,
		long int yPartitionAddr, TInteger64Matrix *yPosAddr, TInteger64Matrix *yRadiiAddr,
		double p, double slack,
		TInteger64Matrix *Pointer) {

	int xres,dim;

	xres=xpos->dimensions[0];
	dim=xpos->dimensions[1];

	// shielding generator setup
	TVarListHandler *xNeighbours;
	xNeighbours=new TVarListHandler();
	xNeighbours->fillFromCSRIndexList(NeighXindices->data, NeighXindptr->data, xres, NeighXindices->dimensions[0]);

	THierarchicalPartition *yPartition=(THierarchicalPartition*) yPartitionAddr;
	double **yPos= (double**) yPosAddr->data;
	double **yRadii = (double**) yRadiiAddr->data;

	TShieldGeneratorTree_PEuclidean *shieldGenerator;
	shieldGenerator= new TShieldGeneratorTree_PEuclidean(dim, yPartition, yPos, yRadii,
			lBottom, lTop,
			xpos->data, xNeighbours, p, slack);

	Pointer->data[0]=(long int) shieldGenerator;

	return 0;
}


int Setup_Shielding_TreeSphere(
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		int lBottom, int lTop,
		TDoubleMatrix *xpos,
		long int yPartitionAddr, TInteger64Matrix *yPosAddr, TInteger64Matrix *yRadiiAddr,
		double p,
		TInteger64Matrix *Pointer) {


	int xres,dim;

	xres=xpos->dimensions[0];
	dim=xpos->dimensions[1];

	// shielding generator setup
	TVarListHandler *xNeighbours;
	xNeighbours=new TVarListHandler();
	xNeighbours->fillFromCSRIndexList(NeighXindices->data, NeighXindptr->data, xres, NeighXindices->dimensions[0]);


	dim=xpos->dimensions[1];
	THierarchicalPartition *yPartition=(THierarchicalPartition*) yPartitionAddr;
	double **yPos= (double**) yPosAddr->data;
	double **yRadii = (double**) yRadiiAddr->data;

	TShieldGeneratorTree_Sphere *shieldGenerator;
	shieldGenerator= new TShieldGeneratorTree_Sphere(dim, yPartition, yPos, yRadii,
			lBottom, lTop,
			xpos->data, xNeighbours, p);

	Pointer->data[0]=(long int) shieldGenerator;

	return 0;
}



int Setup_Shielding_TreeReflector(
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		int lBottom, int lTop,
		TDoubleMatrix *xpos,
		long int yPartitionAddr, TInteger64Matrix *yPosAddr, TInteger64Matrix *yRadiiAddr,
		TInteger64Matrix *Pointer) {


	int xres,dim;

	xres=xpos->dimensions[0];
	dim=xpos->dimensions[1];

	// shielding generator setup
	TVarListHandler *xNeighbours;
	xNeighbours=new TVarListHandler();
	xNeighbours->fillFromCSRIndexList(NeighXindices->data, NeighXindptr->data, xres, NeighXindices->dimensions[0]);


	dim=xpos->dimensions[1];
	THierarchicalPartition *yPartition=(THierarchicalPartition*) yPartitionAddr;
	double **yPos= (double**) yPosAddr->data;
	double **yRadii = (double**) yRadiiAddr->data;

	TShieldGeneratorTree_Reflector *shieldGenerator;
	shieldGenerator= new TShieldGeneratorTree_Reflector(dim, yPartition, yPos, yRadii,
			lBottom, lTop,
			xpos->data, xNeighbours);

	Pointer->data[0]=(long int) shieldGenerator;

	return 0;
}

/* Interaction with shielding generator */


int Shielding_Generate_Request(long int shieldGeneratorAddr, TInteger32Matrix *indices, TInteger32Matrix *indptr,
		TInteger64Matrix *shieldPointer, TInteger32Matrix *shieldSpecs) {

	TVarListHandler *varList,*varListCopy;
	varList=SetupXVars(indices->data, indptr->data, indptr->dimensions[0]-1, indices->dimensions[0]);
	varListCopy=new TVarListHandler(varList);

	TShieldGeneratorBase *shieldGenerator;
	shieldGenerator=(TShieldGeneratorBase*) shieldGeneratorAddr;

	shieldGenerator->generateShield(varList,varListCopy);
	varList->sort();

	shieldPointer->data[0]=(long int) varList;
	shieldSpecs->data[0]=varList->res;
	shieldSpecs->data[1]=varList->total;

	delete varListCopy;
	return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////


int Setup_Shielding_Tree_Benchmark(
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		int lBottom, int lTop,
		TDoubleMatrix *xpos,
		long int yPartitionAddr, TInteger64Matrix *yPosAddr, TInteger64Matrix *yRadiiAddr,
		TInteger64Matrix *Pointer) {

	//int xres,yres,dim,msg;
	int xres,dim;

	xres=xpos->dimensions[0];
	dim=xpos->dimensions[1];

	// shielding generator setup
	TVarListHandler *xNeighbours;
	xNeighbours=new TVarListHandler();
	xNeighbours->fillFromCSRIndexList(NeighXindices->data, NeighXindptr->data, xres, NeighXindices->dimensions[0]);


	THierarchicalPartition *yPartition=(THierarchicalPartition*) yPartitionAddr;
	double **yPos= (double**) yPosAddr->data;
	double **yRadii = (double**) yRadiiAddr->data;

	TShieldGeneratorTree_SqrEuclideanPrototype<TShieldGeneratorTreeBase_Benchmark> *shieldGenerator;
	shieldGenerator= new TShieldGeneratorTree_SqrEuclideanPrototype<TShieldGeneratorTreeBase_Benchmark>(dim, yPartition, yPos, yRadii,
			lBottom, lTop,
			xpos->data, xNeighbours);


	Pointer->data[0]=(long int) shieldGenerator;

	return 0;
}



int Setup_Shielding_TreePEucl_Benchmark(
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		int lBottom, int lTop,
		TDoubleMatrix *xpos,
		long int yPartitionAddr, TInteger64Matrix *yPosAddr, TInteger64Matrix *yRadiiAddr,
		double p, double slack,
		TInteger64Matrix *Pointer) {

	int xres,dim;

	xres=xpos->dimensions[0];
	dim=xpos->dimensions[1];

	// shielding generator setup
	TVarListHandler *xNeighbours;
	xNeighbours=new TVarListHandler();
	xNeighbours->fillFromCSRIndexList(NeighXindices->data, NeighXindptr->data, xres, NeighXindices->dimensions[0]);

	THierarchicalPartition *yPartition=(THierarchicalPartition*) yPartitionAddr;
	double **yPos= (double**) yPosAddr->data;
	double **yRadii = (double**) yRadiiAddr->data;

	TShieldGeneratorTree_PEuclideanPrototype<TShieldGeneratorTreeBase_Benchmark> *shieldGenerator;
	shieldGenerator= new TShieldGeneratorTree_PEuclideanPrototype<TShieldGeneratorTreeBase_Benchmark>(dim, yPartition, yPos, yRadii,
			lBottom, lTop,
			xpos->data, xNeighbours, p, slack);

	Pointer->data[0]=(long int) shieldGenerator;

	return 0;
}


int Setup_Shielding_TreeSphere_Benchmark(
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		int lBottom, int lTop,
		TDoubleMatrix *xpos,
		long int yPartitionAddr, TInteger64Matrix *yPosAddr, TInteger64Matrix *yRadiiAddr,
		double p,
		TInteger64Matrix *Pointer) {


	int xres,dim;

	xres=xpos->dimensions[0];
	dim=xpos->dimensions[1];

	// shielding generator setup
	TVarListHandler *xNeighbours;
	xNeighbours=new TVarListHandler();
	xNeighbours->fillFromCSRIndexList(NeighXindices->data, NeighXindptr->data, xres, NeighXindices->dimensions[0]);


	dim=xpos->dimensions[1];
	THierarchicalPartition *yPartition=(THierarchicalPartition*) yPartitionAddr;
	double **yPos= (double**) yPosAddr->data;
	double **yRadii = (double**) yRadiiAddr->data;

	TShieldGeneratorTree_SpherePrototype<TShieldGeneratorTreeBase_Benchmark> *shieldGenerator;
	shieldGenerator= new TShieldGeneratorTree_SpherePrototype<TShieldGeneratorTreeBase_Benchmark>(dim, yPartition, yPos, yRadii,
			lBottom, lTop,
			xpos->data, xNeighbours, p);

	Pointer->data[0]=(long int) shieldGenerator;

	return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
/* Setup Solver */


int Setup_ShortCutSolver(long int couplingHandlerInterfaceAddr, long int solverInterfaceAddr,
		long int shieldGeneratorAddr, int checkMethod, TInteger64Matrix *SolverPointer) {

	int msg;

	TCouplingHandlerExtBase *couplingHandler;
	TSolverInterface *solverInterface;
	TShieldGeneratorBase *shieldGenerator;

	couplingHandler=(TCouplingHandlerExtBase*) couplingHandlerInterfaceAddr;
	solverInterface=(TSolverInterface*) solverInterfaceAddr;
	shieldGenerator=(TShieldGeneratorBase*) shieldGeneratorAddr;

	TShortCutSolver *ShortCutSolver;
	ShortCutSolver= new TShortCutSolver(couplingHandler,solverInterface,shieldGenerator, checkMethod);

	SolverPointer->data[0]=(long int) ShortCutSolver;


	msg=ShortCutSolver->initialize();
	return msg;
}


//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
/* Interact with Solver */


int ShortCutSolverStep(long int ShortCutSolverAddr, int steps) {
	TShortCutSolver *ShortCutSolver=(TShortCutSolver*) ShortCutSolverAddr;
	int msg;
	msg=ShortCutSolver->step(steps);
	return msg;

}

TShortCutSolverReport ShortCutSolverGetReport(long int ShortCutSolverAddr) {
	TShortCutSolver *ShortCutSolver=(TShortCutSolver*) ShortCutSolverAddr;
	return ShortCutSolver->report;
}


int ShortCutSolverGetSupport_Request(long int ShortCutSolverAddr, TInteger64Matrix *SupportPointer,
		TInteger32Matrix *SupportSpecs) {
	TShortCutSolver *ShortCutSolver=(TShortCutSolver*) ShortCutSolverAddr;
	TVarListHandler *support;
	support=ShortCutSolver->couplingHandler->getSupport();

	SupportPointer->data[0]=(long int) support;
	SupportSpecs->data[0]=support->res;
	SupportSpecs->data[1]=support->total;

	return 0;
}

int ShortCutSolverGetXVars_Request(long int ShortCutSolverAddr, TInteger64Matrix *SupportPointer,
		TInteger32Matrix *SupportSpecs) {
	TShortCutSolver *ShortCutSolver=(TShortCutSolver*) ShortCutSolverAddr;
	TVarListHandler *xVars;
	xVars=ShortCutSolver->couplingHandler->getXVars();

	SupportPointer->data[0]=(long int) xVars;
	SupportSpecs->data[0]=xVars->res;
	SupportSpecs->data[1]=xVars->total;

	return 0;
}

int ShortCutSolverClose(long int ShortCutSolverAddr) {
	TShortCutSolver *ShortCutSolver=(TShortCutSolver*) ShortCutSolverAddr;
	delete ShortCutSolver->couplingHandler;
	delete ShortCutSolver->solverInterface;
	delete ShortCutSolver->shieldGenerator;
	delete ShortCutSolver;
	return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
// Additional Support and Tool Functions

int Tools_Collect_VarList(TInteger32Matrix *indices, TInteger32Matrix *indptr,
		long int varListPointer, int doDelete) {
	TVarListHandler *varList;
	varList=(TVarListHandler*) varListPointer;
	varList->writeToCSRIndexList(indices->data,indptr->data);
	if (doDelete>0) {
		delete varList;
	}
	return 0;
}


int Tools_GetGridNeighbours_Request(TInteger32Matrix *xDims, TInteger64Matrix *NeighbourPointer,
		TInteger32Matrix *NeighbourSpecs, int torusDim) {
	TVarListHandler *neighbours;
	int xres;
	int dim,i;
	dim=xDims->dimensions[0];

	xres=1;
	for(i=0;i<dim;i++) {
		xres=xres*xDims->data[i];
	}

	neighbours=new TVarListHandler();
	neighbours->setupEmpty(xres);

	GridToolsGetNeighbours_Torus(dim, xDims->data, torusDim, neighbours);

	NeighbourPointer->data[0]=(long int) neighbours;
	NeighbourSpecs->data[0]=xres;
	NeighbourSpecs->data[1]=neighbours->total;

	return 0;
}


int Test_VerifyShielding(TDoubleMatrix *c,
		TInteger32Matrix *VLXindices, TInteger32Matrix *VLXindptr,
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		TInteger32Matrix *xMap,
		TInteger64Matrix *missesPointer, TInteger32Matrix *missesSpecs
		) {

	int xres,yres;

	xres=c->dimensions[0];
	yres=c->dimensions[1];

	TVarListHandler *xVars;
	xVars=SetupXVars(VLXindices->data, VLXindptr->data, xres, VLXindices->dimensions[0]);

	// shielding generator setup
	TVarListHandler *xNeighbours;
	xNeighbours=new TVarListHandler();
	xNeighbours->fillFromCSRIndexList(NeighXindices->data, NeighXindptr->data, xres, NeighXindices->dimensions[0]);

	TShieldingVerification shieldingVerification(c->data,xres,yres,xNeighbours);
	TVarListHandler *result;
	result=shieldingVerification.verify(xVars,xMap->data);

	missesPointer->data[0]=(long int) result;
	missesSpecs->data[0]=result->res;
	missesSpecs->data[1]=result->total;

	delete xNeighbours;
	delete xVars;

	return (int) (result->total==0);
}

int Test_VerifyShieldingDuplex(TDoubleMatrix *c,
		TInteger32Matrix *VLXindices, TInteger32Matrix *VLXindptr,
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		TInteger32Matrix *NeighYindices, TInteger32Matrix *NeighYindptr,
		TInteger32Matrix *SPTXindices, TInteger32Matrix *SPTXindptr,
		TInteger64Matrix *missesPointer, TInteger32Matrix *missesSpecs
		) {
	int xres,yres;

	xres=c->dimensions[0];
	yres=c->dimensions[1];

	TVarListHandler *xVars;
	xVars=SetupXVars(VLXindices->data, VLXindptr->data, xres, VLXindices->dimensions[0]);

	// shielding generator setup
	TVarListHandler *xNeighbours,*yNeighbours;
	xNeighbours=new TVarListHandler();
	xNeighbours->fillFromCSRIndexList(NeighXindices->data, NeighXindptr->data, xres, NeighXindices->dimensions[0]);
	yNeighbours=new TVarListHandler();
	yNeighbours->fillFromCSRIndexList(NeighYindices->data, NeighYindptr->data, xres, NeighYindices->dimensions[0]);

	TVarListHandler *xSupport,*ySupport;
	xSupport=SetupXVars(SPTXindices->data, SPTXindptr->data, xres, SPTXindices->dimensions[0]);
	ySupport=new TVarListHandler;
	ySupport->fillViaTranspose(xSupport,yres);

	TShieldingVerificationDuplex shieldingVerification(c->data,xres,yres,xNeighbours,yNeighbours);
	TVarListHandler *result;

	result=shieldingVerification.verify(xVars,xSupport,ySupport);

	missesPointer->data[0]=(long int) result;
	missesSpecs->data[0]=result->res;
	missesSpecs->data[1]=result->total;



	delete xNeighbours;
	delete yNeighbours;
	delete xVars;
	delete xSupport;
	delete ySupport;

	return (int) (result->total==0);
}



int Test_VerifyDualConstraints(TDoubleMatrix *c, TDoubleMatrix *alpha, TDoubleMatrix *beta, double slack,
		TInteger64Matrix *missesPointer, TInteger32Matrix *missesSpecs) {
	int xres,yres,x,y;

	xres=c->dimensions[0];
	yres=c->dimensions[1];

	// varList to capture all violated constraints.
	TVarListHandler *violations;
	violations=new TVarListHandler();
	violations->setupEmpty(xres);


	for(x=0;x<xres;x++) {
		for(y=0;y<yres;y++) {
			if(c->data[x*yres+y]-alpha->data[x]-beta->data[y]<=-slack) {
				// if a violated constraint is found, add to var list
				violations->addToLine(x,y);
			}
		}
	}

	missesPointer->data[0]=(long int) violations;
	missesSpecs->data[0]=violations->res;
	missesSpecs->data[1]=violations->total;

	return (int) (violations->total==0);
}

/*
int Tools_SparseMu_Request(TInteger32Matrix *muSpecs, long int couplingHandlerAddr) {
	TCouplingHandlerSparse *couplingHandler = (TCouplingHandlerSparse*) couplingHandlerAddr;
	muSpecs->data[0]=couplingHandler->total;
	return 0;
}

int Tools_SparseMu_Collect(TDoubleMatrix *mu, long int couplingHandlerAddr, int doDelete) {
	TCouplingHandlerSparse *couplingHandler = (TCouplingHandlerSparse*) couplingHandlerAddr;
	for(int i=0;i<mu->dimensions[0];i++) {
			mu->data[i]=couplingHandler->mu[i];
	}

	if (doDelete>0) {
		delete couplingHandler->xVars;
		delete couplingHandler;
	}
	return 0;
}

*/
