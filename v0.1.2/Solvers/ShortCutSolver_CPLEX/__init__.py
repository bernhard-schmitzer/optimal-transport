import os
import numpy as np
import ctypes as ct
import CClasses
import datetime

_package_directory=os.path.dirname(os.path.abspath(__file__))

_libShortCutSolver = ct.cdll.LoadLibrary(_package_directory+'/Release/libShortCutSolver_CPLEX.so')

CH_SemiDense=0
CH_Sparse=1

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################


def _PreSetupCPLEX(muX, muY, alpha, beta, doMuYnegate, initializeBases):

	# do muY inversion if required
	if doMuYnegate:
		muYParam=np.empty(muY.shape,dtype=np.double)
		muYParam[...]=-muY
	else:
		muYParam=muY

	# allocate dual potential space if required
	if alpha is None:
		alphaParam=np.zeros(muX.shape,dtype=np.double)
	else:
		alphaParam=alpha

	if beta is None:
		betaParam=np.zeros(muY.shape,dtype=np.double)
	else:
		betaParam=beta


	if initializeBases:
		initializeBasesParam=1
	else:
		initializeBasesParam=0

	return (muYParam,alphaParam,betaParam,initializeBasesParam)

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
# Modularized Setup Routines

def Setup_Solver_CPLEX(muX, muY, couplingHandlerPointer, couplingHandlerType, alpha=None, beta=None,
	doMuYnegate=True, initializeBases=False):

	(muYParam,alphaParam,betaParam,initializeBasesParam)=_PreSetupCPLEX(muX, muY, alpha, beta, doMuYnegate, initializeBases)

	pointer=np.zeros((1),dtype=np.int64)
		
	status=_libShortCutSolver.Setup_Solver_CPLEX(\
		CClasses._toCDoubleMatrix(muX),CClasses._toCDoubleMatrix(muYParam),\
		CClasses._toCDoubleMatrix(alphaParam), CClasses._toCDoubleMatrix(betaParam),\
		ct.c_int(initializeBasesParam),\
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


	
def getMethodSetup_SubSolver_CPLEX(couplingHandlerType, initializeBases=True):\
	return lambda muX, muY, couplingHandlerPointer :\
		Setup_Solver_CPLEX(muX, muY, couplingHandlerPointer, couplingHandlerType,\
			doMuYnegate=True, initializeBases=initializeBases)

