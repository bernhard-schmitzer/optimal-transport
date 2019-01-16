#ifndef SINKHORN_H_
#define SINKHORN_H_

#include<iostream>
#include<stdlib.h>
#include<PythonTypes.h>
//#include"TVarListHandler.h"
#include"THierarchicalPartition.h"
#include"TSinkhorn.h"
#include"TMultiVarListHandler.h"


using namespace std;

extern "C" {

long int Setup_CostFunctionProvider_SquaredEuclidean(TInteger64Matrix *posAddr,
		TInteger64Matrix *radiiAddr, TInteger64Matrix *alphaAddr, int posDim, int lBottom, double weight);
long int Setup_CostFunctionProvider_Color_SquaredEuclidean_RGB(TInteger64Matrix *posAddr,
		TInteger64Matrix *radiiAddr, TInteger64Matrix *alphaAddr, int posDim, int lBottom,
		int lTop, double colorWeight, int FR_mode, double FR_kappa);
long int Setup_CostFunctionProvider_Color_SquaredEuclidean_HSV_HS(TInteger64Matrix *posAddr,
		TInteger64Matrix *radiiPosAddr, TInteger64Matrix *radiiValAddr, TInteger64Matrix *radiiHueAddr,
		TInteger64Matrix *alphaAddr, int posDim, int lBottom,
		int liftMode, double colorWeight, int FR_mode, double FR_kappa);
long int Setup_CostFunctionProvider_SquaredEuclideanWF(TInteger64Matrix *posAddr,
		TInteger64Matrix *radiiAddr, TInteger64Matrix *alphaAddr, int posDim, int lBottom, double delta, double cMax);
long int Setup_CostFunctionProvider_SquaredEuclideanBarycenter(TInteger64Matrix *posAddr,
		TInteger64Matrix *radiiAddr, TInteger64Matrix *alphaAddr, TDoubleMatrix *lambda, int posDim, int lBottom);
long int Setup_CostFunctionProvider_Coulomb(TInteger64Matrix *posAddr,
		TInteger64Matrix *radiiAddr, TInteger64Matrix *alphaAddr, int posDim, int lBottom, TDoubleMatrix *charges);
long int Setup_CostFunctionProvider_Reflector_Spherical(TInteger64Matrix *posAddr,
		TInteger64Matrix *radiiAddr, TInteger64Matrix *alphaAddr, int posDim, int lBottom);



long int Setup_CostFunctionProvider_Interpolator(long int coarseAddr, long int fineAddr,
		TInteger64Matrix *partitionAddr, double q, TInteger64Matrix *alphaAddr);

int Check_DualConstraints_Pos(TInteger64Matrix *PartitionAddr,
		int lBottom,
		double slack,
		long int costFunctionProviderAddr,
		TInteger32Matrix *varListSpecs, TInteger64Matrix *varListPointer);


int Refine_VarList_CSR_Pos(TInteger64Matrix *PartitionAddr,
		int lTop,
		long int costFunctionProviderAddr,
		TInteger32Matrix *indices, TInteger32Matrix *indptr,
		TInteger32Matrix *varListSpecs, TInteger64Matrix *varListPointer);

int Refine_VarList_Pos_Pos(TInteger64Matrix *PartitionAddr,
		int lTop,
		long int costFunctionProviderAddr,
		long int kernelAddr,
		TInteger32Matrix *varListSpecs, TInteger64Matrix *varListPointer);


int ReEvaluate_VarList_Pos(TInteger64Matrix *PartitionAddr,
		int lTop,
		long int costFunctionProviderAddr,
		TInteger32Matrix *indices, TInteger32Matrix *indptr,
		TInteger32Matrix *varListSpecs, TInteger64Matrix *varListPointer);


int Iterate(long int kernelAddr, TInteger64Matrix *scalingAddr, TInteger64Matrix *muAddr,
		TInteger32Matrix *res, int n);

int Tools_Collect_PosVarList(TDoubleMatrix *signal, TInteger32Matrix *indices,
		long int kernelPointer, int doDelete);
int Tools_Export_PosVarList(TDoubleMatrix *signal, TInteger32Matrix *indices,
		TInteger64Matrix *pointer);
int Tools_GetSpecs_PosVarList(long int kernelPointer, TInteger32Matrix *specs);
int Tools_Delete_PosVarList(long int kernelPointer);


int Tools_Delete_CostFunctionProvider(long int CFPAddr);
int Tools_Evaluate_CostFunctionProvider(long int CFPAddr, TInteger32Matrix *pos, TDoubleMatrix *c, int layer);


int Tools_GetMarginal(long int kernelAddr, TInteger64Matrix *scalingAddr, TDoubleMatrix *mu,
		TInteger32Matrix *res, int axis);


int Tools_GetDenseCosts_Pos(TInteger64Matrix *PartitionAddr, int layer,
		long int costFunctionProviderAddr,
		TInteger32Matrix *varListSpecs, TInteger64Matrix *varListPointer);


int Kernel_Scale(long int kernelAddr, TInteger64Matrix *muAddr);
int Kernel_Exponentiate(long int kernelAddr, double eps);
}



#endif
