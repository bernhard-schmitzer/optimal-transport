import os
import numpy as np
import ctypes as ct
import CClasses

_package_directory=os.path.dirname(os.path.abspath(__file__))


_libCPLEX = ct.cdll.LoadLibrary(_package_directory+'/Release/libOT_CPLEX.so')

def solveDense(c,muX,muY,alpha=None,beta=None,mu=None,doMuYnegate=True,getDualVariables=False):
	"""Uses the CPLEX network flow solver to solve a dense optimal transport problem between muX and muY, w.r.t. cost function c.
	If mu is given, optimal coupling will be stored in mu. Otherwise, new memory will be allocated.
	If doMuYnegate=True, then a copy of muY is created and negated, to be handed to the solver (as required).
	If False, this has already been done."""
	
	# allocate coupling space if required
	if mu is None:
		muParam=np.zeros(c.shape,dtype=np.double)
	else:
		muParam=mu
	
	# do muY inversion if required
	if doMuYnegate:
		muYParam=np.empty(muY.shape,dtype=np.double)
		muYParam[...]=-muY
	else:
		muYParam=muY
	
	if getDualVariables:
		# if dual variables are requested
	
		# allocate dual potential space if required
		if alpha is None:
			alphaParam=np.zeros(muX.shape,dtype=np.double)
		else:
			alphaParam=alpha

		if beta is None:
			betaParam=np.zeros(muY.shape,dtype=np.double)
		else:
			betaParam=beta

		status=_libCPLEX.solveDense(CClasses._toCDoubleMatrix(c), CClasses._toCDoubleMatrix(muX),\
			CClasses._toCDoubleMatrix(muYParam), CClasses._toCDoubleMatrix(muParam),True,\
			CClasses._toCDoubleMatrix(alphaParam),CClasses._toCDoubleMatrix(betaParam)\
			)
	
		result={"status" : status, "mu" : muParam, "alpha" : alphaParam, "beta" : betaParam}
		return result

	else:
		status=_libCPLEX.solveDense(CClasses._toCDoubleMatrix(c), CClasses._toCDoubleMatrix(muX),\
			CClasses._toCDoubleMatrix(muYParam), CClasses._toCDoubleMatrix(muParam),False,0,0\
			)
	
		result={"status" : status, "mu" : muParam}
		return result

def solveSemiDense(c,muX,muY,indices,sep,alpha=None,beta=None,mu=None,doMuYnegate=True,getDualVariables=False):
	"""Uses the Hungarian method to solve a dense optimal transport problem between muX and muY, w.r.t. cost function c.
	Uses semi sparse coupling handler. (indices,sep) specify location of active variables in scipy.sparse.csr_matrix format.
	If mu is given, optimal coupling will be stored in mu. Otherwise, new memory will be allocated.
	If doMuYnegate=True, then a copy of muY is created and negated, to be handed to the solver (as required).
	If False, this has already been done."""

	if mu is None:
		muParam=np.zeros(c.shape,dtype=np.double)
	else:
		muParam=mu
	
	if doMuYnegate:
		muYParam=np.empty(muY.shape,dtype=np.double)
		muYParam[...]=-muY
	else:
		muYParam=muY

	if getDualVariables:
		# if dual variables are requested
	
		# allocate dual potential space if required
		if alpha is None:
			alphaParam=np.zeros(muX.shape,dtype=np.double)
		else:
			alphaParam=alpha

		if beta is None:
			betaParam=np.zeros(muY.shape,dtype=np.double)
		else:
			betaParam=beta

		status=_libCPLEX.solveSemiDense(CClasses._toCDoubleMatrix(c), CClasses._toCDoubleMatrix(muX),
			CClasses._toCDoubleMatrix(muYParam),\
			CClasses._toCDoubleMatrix(muParam),\
			True,CClasses._toCDoubleMatrix(alphaParam),CClasses._toCDoubleMatrix(betaParam),\
			CClasses._toCInteger32Matrix(indices), CClasses._toCInteger32Matrix(sep))
	
	
		result={"status" : status, "mu" : muParam, "alpha" : alphaParam, "beta" : betaParam}
		return result
	else:		
		status=_libCPLEX.solveSemiDense(CClasses._toCDoubleMatrix(c), CClasses._toCDoubleMatrix(muX),
			CClasses._toCDoubleMatrix(muYParam),\
			CClasses._toCDoubleMatrix(muParam),\
			False,0,0,
			CClasses._toCInteger32Matrix(indices), CClasses._toCInteger32Matrix(sep))
	
	
		result={"status" : status, "mu" : muParam}
		return result


def solveSparse(muX,muY,data,indices,sep,alpha=None,beta=None,doMuYnegate=True,getDualVariables=False):

	if doMuYnegate:
		muYParam=np.empty(muY.shape,dtype=np.double)
		muYParam[...]=-muY
	else:
		muYParam=muY

	if getDualVariables:
		# if dual variables are requested
	
		# allocate dual potential space if required
		if alpha is None:
			alphaParam=np.zeros(muX.shape,dtype=np.double)
		else:
			alphaParam=alpha

		if beta is None:
			betaParam=np.zeros(muY.shape,dtype=np.double)
		else:
			betaParam=beta
		
		dualParamC=True
		alphaParamC=CClasses._toCDoubleMatrix(alphaParam)
		betaParamC=CClasses._toCDoubleMatrix(betaParam)
	else:
		dualParamC=False
		alphaParamC=0
		betaParamC=0


	muSpecs=np.zeros((1),dtype=np.int32)
	muPointer=np.zeros((1),dtype=np.int64)

	status=_libCPLEX.solveSparse(\
		ct.c_int(muX.shape[0]), ct.c_int(muY.shape[0]),\
		CClasses._toCDoubleMatrix(data), CClasses._toCDoubleMatrix(muX), CClasses._toCDoubleMatrix(muYParam),\
		dualParamC,alphaParamC,betaParamC,\
		CClasses._toCInteger32Matrix(indices), CClasses._toCInteger32Matrix(sep),\
		CClasses._toCInteger64Matrix(muSpecs), CClasses._toCInteger64Matrix(muPointer)\
		)
	
	mu=np.zeros((muSpecs[0]),dtype=np.double)
	_libCPLEX.solveSparse_getMu(CClasses._toCDoubleMatrix(mu),ct.c_long(muPointer[0]))
	
	if getDualVariables:
		result={"status" : status, "alpha" : alphaParam, "beta" : betaParam, "mu" : mu}
	else:
		result={"status" : status, "mu" : mu}

	return result

def solveSparse_fullC(c,muX,muY,indices,sep,alpha=None,beta=None,doMuYnegate=True,getDualVariables=False):

	if doMuYnegate:
		muYParam=np.empty(muY.shape,dtype=np.double)
		muYParam[...]=-muY
	else:
		muYParam=muY

	if getDualVariables:
		# if dual variables are requested
	
		# allocate dual potential space if required
		if alpha is None:
			alphaParam=np.zeros(muX.shape,dtype=np.double)
		else:
			alphaParam=alpha

		if beta is None:
			betaParam=np.zeros(muY.shape,dtype=np.double)
		else:
			betaParam=beta
		
		dualParamC=True
		alphaParamC=CClasses._toCDoubleMatrix(alphaParam)
		betaParamC=CClasses._toCDoubleMatrix(betaParam)
	else:
		dualParamC=False
		alphaParamC=0
		betaParamC=0


	muSpecs=np.zeros((1),dtype=np.int32)
	muPointer=np.zeros((1),dtype=np.int64)

	status=_libCPLEX.solveSparse_fullC(\
		CClasses._toCDoubleMatrix(c),\
		CClasses._toCDoubleMatrix(muX), CClasses._toCDoubleMatrix(muYParam),\
		dualParamC,alphaParamC,betaParamC,\
		CClasses._toCInteger32Matrix(indices), CClasses._toCInteger32Matrix(sep),\
		CClasses._toCInteger64Matrix(muSpecs), CClasses._toCInteger64Matrix(muPointer)\
		)
	
	mu=np.zeros((muSpecs[0]),dtype=np.double)
	_libCPLEX.solveSparse_getMu(CClasses._toCDoubleMatrix(mu),ct.c_long(muPointer[0]))
	
	if getDualVariables:
		result={"status" : status, "alpha" : alphaParam, "beta" : betaParam, "mu" : mu}
	else:
		result={"status" : status, "mu" : mu}

	return result
