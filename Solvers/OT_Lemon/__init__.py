import os
import numpy as np
import ctypes as ct
import CClasses

_package_directory=os.path.dirname(os.path.abspath(__file__))

_libSolverLemon = ct.cdll.LoadLibrary(_package_directory+'/Release/libOT_Lemon.so')




##########################################################################################################################
# Simple Solvers


def solveDense(c,muX,muY,cScale,algorithm=0,alpha=None,beta=None,mu=None):

	# allocate coupling space if required
	if mu is None:
		muParam=np.zeros(c.shape,dtype=np.double)
	else:
		muParam=mu
	
	# allocate dual potential space if required
	if alpha is None:
		alphaParam=np.zeros(muX.shape,dtype=np.double)
	else:
		alphaParam=alpha

	if beta is None:
		betaParam=np.zeros(muY.shape,dtype=np.double)
	else:
		betaParam=beta
	
	
	status=_libSolverLemon.solveDense(CClasses._toCDoubleMatrix(c), CClasses._toCDoubleMatrix(muX),\
		CClasses._toCDoubleMatrix(muY), CClasses._toCDoubleMatrix(muParam),True,\
		CClasses._toCDoubleMatrix(alphaParam),CClasses._toCDoubleMatrix(betaParam),\
		ct.c_double(cScale), ct.c_int(algorithm)\
		)
	
	result={"status" : status, "mu" : muParam, "alpha" : alphaParam, "beta" : betaParam}
	return result


def solveSemiDense(c,muX,muY,indices,sep,cScale,algorithm=0,alpha=None,beta=None,mu=None):

	# allocate coupling space if required
	if mu is None:
		muParam=np.zeros(c.shape,dtype=np.double)
	else:
		muParam=mu
	
	# allocate dual potential space if required
	if alpha is None:
		alphaParam=np.zeros(muX.shape,dtype=np.double)
	else:
		alphaParam=alpha

	if beta is None:
		betaParam=np.zeros(muY.shape,dtype=np.double)
	else:
		betaParam=beta
	
	status=_libSolverLemon.solveDense(CClasses._toCDoubleMatrix(c), CClasses._toCDoubleMatrix(muX),\
		CClasses._toCDoubleMatrix(muY), CClasses._toCDoubleMatrix(muParam),True,\
		CClasses._toCDoubleMatrix(alphaParam),CClasses._toCDoubleMatrix(betaParam),\
		CClasses._toCInteger32Matrix(indices), CClasses._toCInteger32Matrix(sep),\
		ct.c_double(cScale), ct.c_int(algorithm)\
		)
	
	result={"status" : status, "mu" : muParam, "alpha" : alphaParam, "beta" : betaParam}
	return result


