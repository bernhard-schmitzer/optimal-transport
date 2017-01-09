import os
import numpy as np
import ctypes as ct
import CClasses
import datetime

_package_directory=os.path.dirname(os.path.abspath(__file__))

_libShortCutSolver = ct.cdll.LoadLibrary(_package_directory+'/Release/libShortCutSolver_Lemon.so')

CH_SemiDense=0
CH_Sparse=1

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################


def _PreSetupLemon(muX, muY, alpha, beta, dualOffset):

	# allocate dual potential space if required
	if alpha is None:
		alphaParam=np.zeros(muX.shape,dtype=np.double)
	else:
		alphaParam=alpha

	if beta is None:
		betaParam=np.zeros(muY.shape,dtype=np.double)
	else:
		betaParam=beta


	if dualOffset:
		dualOffsetParam=1
	else:
		dualOffsetParam=0

	return (alphaParam,betaParam,dualOffsetParam)

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
# Modularized Setup Routines

#int Setup_Solver_Lemon(TDoubleMatrix *muX, TDoubleMatrix *muY,
#		TDoubleMatrix *alpha, TDoubleMatrix *beta,
#		int dualOffset, double cScale, int algorithm,
#		long int couplingHandlerAddr,
#		int couplingHandlerType,
#		TInteger64Matrix *Pointer);


def Setup_Solver_Lemon(muX, muY, couplingHandlerPointer, couplingHandlerType, cScale, algorithm=0, alpha=None, beta=None,
	dualOffset=False):

	(alphaParam,betaParam,dualOffsetParam)=_PreSetupLemon(muX, muY, alpha, beta, dualOffset)

	pointer=np.zeros((1),dtype=np.int64)
		
	status=_libShortCutSolver.Setup_Solver_Lemon(\
		CClasses._toCDoubleMatrix(muX),CClasses._toCDoubleMatrix(muY),\
		CClasses._toCDoubleMatrix(alphaParam), CClasses._toCDoubleMatrix(betaParam),\
		ct.c_int(dualOffsetParam), ct.c_double(cScale), ct.c_int(algorithm),\
		ct.c_long(couplingHandlerPointer),\
		ct.c_int(couplingHandlerType),\
		CClasses._toCInteger64Matrix(pointer)\
		)
	
	result={"status" : status, "pointer" : pointer[0], "alpha" : alphaParam, "beta" : betaParam}
	return result


########################################################################################################################
########################################################################################################################
########################################################################################################################
# Prepare setup methods


	
def getMethodSetup_SubSolver_Lemon(couplingHandlerType, cScale, algorithm=0, dualOffset=True):\
	return lambda muX, muY, couplingHandlerPointer :\
		Setup_Solver_Lemon(muX, muY, couplingHandlerPointer, couplingHandlerType, cScale=cScale,\
			algorithm=algorithm, dualOffset=dualOffset)

