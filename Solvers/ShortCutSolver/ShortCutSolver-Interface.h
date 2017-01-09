#ifndef ShortCutSolver_Interface_H
#define ShortCutSolver_Interface_H

#include<stdlib.h>
#include<PythonTypes.h>
#include"TCouplingHandler.h"
#include"TCostFunctionProvider-Dynamic.h"
#include"TVarListHandler.h"
#include"THierarchicalPartition.h"
#include"ShieldGenerator/TShieldGenerator.h"
#include"TShortCutSolver.h"
#include"ShortCutSolver-Tools.h"



using namespace std;


extern "C" {

//////////////////////////////////////////////////////////////////////////////////////////////
/* Cost Function Providers */

int Setup_CostFunctionProvider_SqrEuclidean(
		TDoubleMatrix *xPos, TDoubleMatrix *yPos,
		TInteger64Matrix *Pointer);

int Setup_CostFunctionProvider_Torus(
		TDoubleMatrix *xPos, TDoubleMatrix *yPos, TDoubleMatrix *radius, int torusDim,
		TInteger64Matrix *Pointer);


int CostFunctionProvider_Evaluate(long int costFunctionProviderAddr,
		 	TInteger32Matrix *VLXindices, TInteger32Matrix *VLXindptr, TDoubleMatrix *c);

int CostFunctionProvider_GetDenseCost(long int costFunctionProviderAddr,
 		 	TDoubleMatrix *c);

int CostFunctionProvider_GetRes(long int costFunctionProviderAddr,
		 	TInteger32Matrix *res);

int CostFunctionProvider_Delete(long int costFunctionProviderAddr);

//////////////////////////////////////////////////////////////////////////////////////////////
/* Coupling Handlers */

int Setup_CouplingHandler_SemiDense(TDoubleMatrix *c, TDoubleMatrix *mu,
		TInteger32Matrix *VLXindices, TInteger32Matrix *VLXindptr,
		TInteger64Matrix *Pointer);

int Setup_CouplingHandler_Sparse_fullC(TDoubleMatrix *c,
		TInteger32Matrix *VLXindices, TInteger32Matrix *VLXindptr,
		TInteger64Matrix *Pointer);

int Setup_CouplingHandler_Sparse_dynamicC(
		TInteger32Matrix *VLXindices, TInteger32Matrix *VLXindptr,
		int xres, int yres, long int costFunctionProviderAddr,
		TInteger64Matrix *Pointer);


/* interaction with sparse coupling handler*/

int CouplingHandler_Sparse_GetMu_Request(
		long int couplingHandlerAddr, TInteger32Matrix *MuSpecs);

int CouplingHandler_Sparse_GetMu_Collect(
		long int couplingHandlerAddr, TDoubleMatrix *data, TInteger32Matrix *indices, TInteger32Matrix *indptr);

int CouplingHandler_Sparse_GetSupport_Request(long int couplingHandlerAddr, TInteger64Matrix *SupportPointer,
		TInteger32Matrix *SupportSpecs);


int CouplingHandler_Sparse_GetSupport_Collect(
		long int supportAddr,
		TDoubleMatrix *data, TInteger32Matrix *indices, TInteger32Matrix *indptr);


int CouplingHandler_Sparse_GetCost(
		long int couplingHandlerAddr,
		TDoubleMatrix *data,
		int total);

	//////////////////////////////////////////////////////////////////////////////////////////////
/* Shielding Methods */



int Setup_Shielding_Grid(TInteger32Matrix *xDims, TInteger32Matrix *yDims, TInteger64Matrix *Pointer);

int Setup_Shielding_Padding(TInteger32Matrix *xDims, TInteger32Matrix *yDims, TInteger64Matrix *Pointer, int width);


int Setup_Shielding_Tree(
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		int lBottom, int lTop,
		TDoubleMatrix *xpos,
		long int yPartitionAddr, TInteger64Matrix *yPosAddr, TInteger64Matrix *yRadiiAddr,
		TInteger64Matrix *Pointer);

int Setup_Shielding_Tree_Torus(
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		int lBottom, int lTop,
		TDoubleMatrix *xpos,
		long int yPartitionAddr, TInteger64Matrix *yPosAddr, TInteger64Matrix *yRadiiAddr,
		TInteger64Matrix *yTorusRadiiAddr,
		TDoubleMatrix *torusRadii, int torusDim,
		TInteger64Matrix *Pointer);


int Setup_Shielding_TreeNoise(
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		int lBottom, int lTop,
		TDoubleMatrix *xpos,
		long int yPartitionAddr, TInteger64Matrix *yPosAddr, TInteger64Matrix *yRadiiAddr,
		TInteger64Matrix *c, double eta, double lambda,
		TInteger64Matrix *Pointer);

int Setup_Shielding_TreePEucl(
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		int lBottom, int lTop,
		TDoubleMatrix *xpos,
		long int yPartitionAddr, TInteger64Matrix *yPosAddr, TInteger64Matrix *yRadiiAddr,
		double p, double slack,
		TInteger64Matrix *Pointer);

int Setup_Shielding_TreeSphere(
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		int lBottom, int lTop,
		TDoubleMatrix *xpos,
		long int yPartitionAddr, TInteger64Matrix *yPosAddr, TInteger64Matrix *yRadiiAddr,
		double p,
		TInteger64Matrix *Pointer);


int Setup_Shielding_TreeReflector(
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		int lBottom, int lTop,
		TDoubleMatrix *xpos,
		long int yPartitionAddr, TInteger64Matrix *yPosAddr, TInteger64Matrix *yRadiiAddr,
		TInteger64Matrix *Pointer);


/* interaction with shielding generators */


int Shielding_Generate_Request(long int shieldGeneratorAddr, TInteger32Matrix *indices, TInteger32Matrix *indptr,
		TInteger64Matrix *shieldPointer, TInteger32Matrix *shieldSpecs);

//////////////////////////////////////////////////////////////////////////////////////////////
/* Shielding Benchmark */


int Setup_Shielding_Tree_Benchmark(
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		int lBottom, int lTop,
		TDoubleMatrix *xpos,
		long int yPartitionAddr, TInteger64Matrix *yPosAddr, TInteger64Matrix *yRadiiAddr,
		TInteger64Matrix *Pointer);


int Setup_Shielding_TreePEucl_Benchmark(
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		int lBottom, int lTop,
		TDoubleMatrix *xpos,
		long int yPartitionAddr, TInteger64Matrix *yPosAddr, TInteger64Matrix *yRadiiAddr,
		double p, double slack,
		TInteger64Matrix *Pointer);

int Setup_Shielding_TreeSphere_Benchmark(
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		int lBottom, int lTop,
		TDoubleMatrix *xpos,
		long int yPartitionAddr, TInteger64Matrix *yPosAddr, TInteger64Matrix *yRadiiAddr,
		double p,
		TInteger64Matrix *Pointer);

//////////////////////////////////////////////////////////////////////////////////////////////
/* Setup Solver */


int Setup_ShortCutSolver(long int couplingHandlerInterfaceAddr, long int solverInterfaceAddr,
		long int shieldGeneratorAddr,int checkMethod, TInteger64Matrix *SolverPointer);

//////////////////////////////////////////////////////////////////////////////////////////////
/* Interact with Solver */


int ShortCutSolverStep(long int ShortCutSolverAddr, int steps);

TShortCutSolverReport ShortCutSolverGetReport(long int ShortCutSolverAddr);

int ShortCutSolverGetSupport_Request(long int ShortCutSolverAddr, TInteger64Matrix *SupportPointer,
		TInteger32Matrix *SupportSpecs);

int ShortCutSolverGetXVars_Request(long int ShortCutSolverAddr, TInteger64Matrix *SupportPointer,
		TInteger32Matrix *SupportSpecs);

int ShortCutSolverClose(long int ShortCutSolverAddr);



////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
// Additional Support and Tool Functions

int Tools_Collect_VarList(TInteger32Matrix *indices, TInteger32Matrix *indptr,
		long int varListPointer, int doDelete);

int Tools_GetGridNeighbours_Request(TInteger32Matrix *xDims, TInteger64Matrix *NeighbourPointer,
		TInteger32Matrix *NeighbourSpecs, int torusDim);

int Test_VerifyShielding(TDoubleMatrix *c,
		TInteger32Matrix *VLXindices, TInteger32Matrix *VLXindptr,
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		TInteger32Matrix *xMap,
		TInteger64Matrix *missesPointer, TInteger32Matrix *missesSpecs
		);

int Test_VerifyShieldingDuplex(TDoubleMatrix *c,
		TInteger32Matrix *VLXindices, TInteger32Matrix *VLXindptr,
		TInteger32Matrix *NeighXindices, TInteger32Matrix *NeighXindptr,
		TInteger32Matrix *NeighYindices, TInteger32Matrix *NeighYindptr,
		TInteger32Matrix *SPTXindices, TInteger32Matrix *SPTXindptr,
		TInteger64Matrix *missesPointer, TInteger32Matrix *missesSpecs
		);


int Test_VerifyDualConstraints(TDoubleMatrix *c, TDoubleMatrix *alpha, TDoubleMatrix *beta, double slack,
		TInteger64Matrix *missesPointer, TInteger32Matrix *missesSpecs);
}

#endif

