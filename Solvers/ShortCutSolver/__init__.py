import os
import numpy as np
import ctypes as ct
import CClasses
import datetime

_package_directory=os.path.dirname(os.path.abspath(__file__))

_libShortCutSolver = ct.cdll.LoadLibrary(_package_directory+'/Release/libShortCutSolver.so')

class TShortCutSolverReport(ct.Structure):
	_fields_=[("steps",ct.c_int),("objective",ct.c_double),("solved",ct.c_int),("t_solving",ct.c_int),("t_shielding",ct.c_int),\
		("n_variables",ct.c_int),("n_shielding_queries",ct.c_int)]

_libShortCutSolver.ShortCutSolverGetReport.restype=TShortCutSolverReport

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
# Setup Final ShortCutSolver


def Setup_ShortCutSolver(pointerCouplingHandler, pointerSolver, pointerShieldGenerator, checkMethod=1):
	pointer=np.zeros((1),dtype=np.int64)
		
	status=_libShortCutSolver.Setup_ShortCutSolver(
		ct.c_long(pointerCouplingHandler),\
		ct.c_long(pointerSolver),\
		ct.c_long(pointerShieldGenerator),\
		ct.c_int(checkMethod),\
		CClasses._toCInteger64Matrix(pointer)\
		)
	result={"status" : status, "pointer" : pointer[0]}
	return result


#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
# Modularized Setup Cost Function Provider
def Setup_CostFunctionProvider_SqrEuclidean(xPos, yPos):

	pointer=np.zeros((1),dtype=np.int64)
		
	status=_libShortCutSolver.Setup_CostFunctionProvider_SqrEuclidean(\
		CClasses._toCDoubleMatrix(xPos), CClasses._toCDoubleMatrix(yPos), \
		CClasses._toCInteger64Matrix(pointer)\
		)
	return pointer[0]


def Setup_CostFunctionProvider_Torus(xPos, yPos, radius, torusDim):

	pointer=np.zeros((1),dtype=np.int64)
		
	status=_libShortCutSolver.Setup_CostFunctionProvider_Torus(\
		CClasses._toCDoubleMatrix(xPos), CClasses._toCDoubleMatrix(yPos), \
		CClasses._toCDoubleMatrix(radius),\
		ct.c_int(torusDim),\
		CClasses._toCInteger64Matrix(pointer)\
		)
	return pointer[0]


# Interaction with cost function provider #

def CostFunctionProvider_Evaluate(pointer, indices, indptr):
	c=np.zeros(indices.shape,dtype=np.double)
	status=_libShortCutSolver.CostFunctionProvider_Evaluate(ct.c_long(pointer),\
		 	CClasses._toCInteger32Matrix(indices),\
		 	CClasses._toCInteger32Matrix(indptr),\
		 	CClasses._toCDoubleMatrix(c))
	return c


def CostFunctionProvider_GetDenseCost(pointer):
	res=CostFunctionProvider_GetRes(pointer)
	c=np.zeros(res,dtype=np.double)
	status=_libShortCutSolver.CostFunctionProvider_GetDenseCost(ct.c_long(pointer),\
		 	CClasses._toCDoubleMatrix(c))
	return c

def CostFunctionProvider_GetRes(pointer):
	res=np.zeros((2),dtype=np.int32)
	status=_libShortCutSolver.CostFunctionProvider_GetRes(ct.c_long(pointer),\
		 	CClasses._toCInteger32Matrix(res))
	return res
	
def CostFunctionProvider_Delete(pointer):
	status=_libShortCutSolver.CostFunctionProvider_Delete(ct.c_long(pointer))


#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
# Modularized Setup Coupling Handlers

def Setup_CouplingHandler_SemiDense(c, VLXindices, VLXindptr, mu=None):

	# allocate coupling space if required
	if mu is None:
		muParam=np.zeros(c.shape,dtype=np.double)
	else:
		muParam=mu

	pointer=np.zeros((1),dtype=np.int64)
		
	status=_libShortCutSolver.Setup_CouplingHandler_SemiDense(CClasses._toCDoubleMatrix(c),\
		CClasses._toCDoubleMatrix(muParam),\
		CClasses._toCInteger32Matrix(VLXindices),CClasses._toCInteger32Matrix(VLXindptr),\
		CClasses._toCInteger64Matrix(pointer)\
		)

	result={"status" : status, "pointer" : pointer[0], "mu" : muParam}
	return result

def Setup_CouplingHandler_Sparse_fullC(c, VLXindices, VLXindptr):

	pointer=np.zeros((1),dtype=np.int64)
		
	status=_libShortCutSolver.Setup_CouplingHandler_Sparse_fullC(CClasses._toCDoubleMatrix(c),\
		CClasses._toCInteger32Matrix(VLXindices),CClasses._toCInteger32Matrix(VLXindptr),\
		CClasses._toCInteger64Matrix(pointer)\
		)
	result={"status" : status, "pointer" : pointer[0]}
	return result

def Setup_CouplingHandler_Sparse_dynamicC(VLXindices, VLXindptr,\
	xres, yres, CFPpointer):

	pointer=np.zeros((1),dtype=np.int64)
		
	status=_libShortCutSolver.Setup_CouplingHandler_Sparse_dynamicC(\
		CClasses._toCInteger32Matrix(VLXindices),CClasses._toCInteger32Matrix(VLXindptr),\
		ct.c_int(xres), ct.c_int(yres), ct.c_long(CFPpointer),\
		CClasses._toCInteger64Matrix(pointer)\
		)
	result={"status" : status, "pointer" : pointer[0]}
	return result

# Interaction with sparse coupling handler

def CouplingHandler_Sparse_GetMu(pointer):

	specs=np.zeros((2),dtype=np.int32)

	_libShortCutSolver.CouplingHandler_Sparse_GetMu_Request(\
			ct.c_long(pointer), CClasses._toCInteger32Matrix(specs)\
			)

	data=np.zeros((specs[1]),dtype=np.double)
	indices=np.zeros((specs[1]),dtype=np.int32)
	indptr=np.zeros((specs[0]+1),dtype=np.int32)

	_libShortCutSolver.CouplingHandler_Sparse_GetMu_Collect(\
			ct.c_long(pointer), CClasses._toCDoubleMatrix(data),\
			CClasses._toCInteger32Matrix(indices), CClasses._toCInteger32Matrix(indptr)\
			)
			
	return (data,indices,indptr)
	

def CouplingHandler_Sparse_GetSupport(pointer):

	specs=np.zeros((2),dtype=np.int32)
	supportPointer=np.zeros((1),dtype=np.int64)

	_libShortCutSolver.CouplingHandler_Sparse_GetSupport_Request(\
			ct.c_long(pointer), CClasses._toCInteger64Matrix(supportPointer), CClasses._toCInteger32Matrix(specs)\
			)

	data=np.zeros((specs[1]),dtype=np.double)
	indices=np.zeros((specs[1]),dtype=np.int32)
	indptr=np.zeros((specs[0]+1),dtype=np.int32)

	_libShortCutSolver.CouplingHandler_Sparse_GetSupport_Collect(\
			ct.c_long(supportPointer[0]), CClasses._toCDoubleMatrix(data),\
			CClasses._toCInteger32Matrix(indices), CClasses._toCInteger32Matrix(indptr)\
			)
			
	return (data,indices,indptr)


def CouplingHandler_Sparse_GetCost(pointer, total):

	data=np.zeros((total),dtype=np.double)

	_libShortCutSolver.CouplingHandler_Sparse_GetCost(\
			ct.c_long(pointer), \
			CClasses._toCDoubleMatrix(data),\
			ct.c_int(total)\
			)
			
	return data

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
# Modularized Setup Shielding Methods

def Setup_Shielding_Grid(xDims, yDims):

	pointer=np.zeros((1),dtype=np.int64)
		
	status=_libShortCutSolver.Setup_Shielding_Grid(
		CClasses._toCInteger32Matrix(xDims),CClasses._toCInteger32Matrix(yDims),\
		CClasses._toCInteger64Matrix(pointer)\
		)
	result={"status" : status, "pointer" : pointer[0]}
	return result


def Setup_Shielding_Padding(xDims, yDims, width):

	pointer=np.zeros((1),dtype=np.int64)
		
	status=_libShortCutSolver.Setup_Shielding_Padding(
		CClasses._toCInteger32Matrix(xDims),CClasses._toCInteger32Matrix(yDims),\
		CClasses._toCInteger64Matrix(pointer),\
		ct.c_int(width)
		)
	result={"status" : status, "pointer" : pointer[0]}
	return result


def Setup_Shielding_Tree(NeighXindices,NeighXindptr,\
	lBottom, lTop, xpos, pointerY, pointerYpos, pointerYradii):

	pointer=np.zeros((1),dtype=np.int64)
	status=_libShortCutSolver.Setup_Shielding_Tree(\
		CClasses._toCInteger32Matrix(NeighXindices),CClasses._toCInteger32Matrix(NeighXindptr),\
		ct.c_int(lBottom), ct.c_int(lTop),
		CClasses._toCDoubleMatrix(xpos),
		ct.c_long(pointerY), CClasses._toCInteger64Matrix(pointerYpos), CClasses._toCInteger64Matrix(pointerYradii),
		CClasses._toCInteger64Matrix(pointer)\
		)

	result={"status" : status, "pointer" : pointer[0]}
	return result

def Setup_Shielding_Tree_Torus(NeighXindices,NeighXindptr,\
	lBottom, lTop, xpos, pointerY, pointerYpos, pointerYradii, pointerYTorusRadii, torusRadii, torusDim):

	pointer=np.zeros((1),dtype=np.int64)
	status=_libShortCutSolver.Setup_Shielding_Tree_Torus(\
		CClasses._toCInteger32Matrix(NeighXindices),CClasses._toCInteger32Matrix(NeighXindptr),\
		ct.c_int(lBottom), ct.c_int(lTop),
		CClasses._toCDoubleMatrix(xpos),
		ct.c_long(pointerY), CClasses._toCInteger64Matrix(pointerYpos), CClasses._toCInteger64Matrix(pointerYradii),\
		CClasses._toCInteger64Matrix(pointerYTorusRadii),\
		CClasses._toCDoubleMatrix(torusRadii), ct.c_int(torusDim),\
		CClasses._toCInteger64Matrix(pointer)\
		)

	result={"status" : status, "pointer" : pointer[0]}
	return result

def Setup_Shielding_TreeNoise(NeighXindices,NeighXindptr,\
	lBottom, lTop, xpos, pointerY, pointerYpos, pointerYradii,\
	pointerCList, c_eta, c_lambda):

	pointer=np.zeros((1),dtype=np.int64)
	status=_libShortCutSolver.Setup_Shielding_TreeNoise(\
		CClasses._toCInteger32Matrix(NeighXindices),CClasses._toCInteger32Matrix(NeighXindptr),\
		ct.c_int(lBottom), ct.c_int(lTop),\
		CClasses._toCDoubleMatrix(xpos),\
		ct.c_long(pointerY), CClasses._toCInteger64Matrix(pointerYpos), CClasses._toCInteger64Matrix(pointerYradii),\
		CClasses._toCInteger64Matrix(pointerCList), ct.c_double(c_eta), ct.c_double(c_lambda),\
		CClasses._toCInteger64Matrix(pointer)\
		)

	result={"status" : status, "pointer" : pointer[0]}
	return result

def Setup_Shielding_TreePEucl(NeighXindices,NeighXindptr,\
	lBottom, lTop, xpos, pointerY, pointerYpos, pointerYradii, p, slack=0.):

	pointer=np.zeros((1),dtype=np.int64)
	status=_libShortCutSolver.Setup_Shielding_TreePEucl(\
		CClasses._toCInteger32Matrix(NeighXindices),CClasses._toCInteger32Matrix(NeighXindptr),\
		ct.c_int(lBottom), ct.c_int(lTop),\
		CClasses._toCDoubleMatrix(xpos),\
		ct.c_long(pointerY), CClasses._toCInteger64Matrix(pointerYpos), CClasses._toCInteger64Matrix(pointerYradii),\
		ct.c_double(p), ct.c_double(slack),\
		CClasses._toCInteger64Matrix(pointer)\
		)

	result={"status" : status, "pointer" : pointer[0]}
	return result

def Setup_Shielding_TreeSphere(NeighXindices, NeighXindptr,\
	lBottom, lTop, xpos, pointerY, pointerYpos, pointerYradii, p):

	pointer=np.zeros((1),dtype=np.int64)
	status=_libShortCutSolver.Setup_Shielding_TreeSphere(\
		CClasses._toCInteger32Matrix(NeighXindices),CClasses._toCInteger32Matrix(NeighXindptr),\
		ct.c_int(lBottom), ct.c_int(lTop),\
		CClasses._toCDoubleMatrix(xpos),\
		ct.c_long(pointerY), CClasses._toCInteger64Matrix(pointerYpos), CClasses._toCInteger64Matrix(pointerYradii),\
		ct.c_double(p),\
		CClasses._toCInteger64Matrix(pointer)\
		)

	result={"status" : status, "pointer" : pointer[0]}
	return result

def Setup_Shielding_TreeReflector(NeighXindices, NeighXindptr,\
	lBottom, lTop, xpos, pointerY, pointerYpos, pointerYradii):

	pointer=np.zeros((1),dtype=np.int64)
	status=_libShortCutSolver.Setup_Shielding_TreeReflector(\
		CClasses._toCInteger32Matrix(NeighXindices),CClasses._toCInteger32Matrix(NeighXindptr),\
		ct.c_int(lBottom), ct.c_int(lTop),\
		CClasses._toCDoubleMatrix(xpos),\
		ct.c_long(pointerY), CClasses._toCInteger64Matrix(pointerYpos), CClasses._toCInteger64Matrix(pointerYradii),\
		CClasses._toCInteger64Matrix(pointer)\
		)

	result={"status" : status, "pointer" : pointer[0]}
	return result


#####################################################################################################################################
# Interaction with Shield Generator

def Shielding_Generate(pointerShielding, varList):

	specs=np.zeros((2),dtype=np.int32)
	pointerVarList=np.zeros((1),dtype=np.int64)

	_libShortCutSolver.Shielding_Generate_Request(\
			ct.c_long(pointerShielding),
			CClasses._toCInteger32Matrix(varList[0]), CClasses._toCInteger32Matrix(varList[1]), \
			CClasses._toCInteger64Matrix(pointerVarList), CClasses._toCInteger32Matrix(specs)\
			)

	indices=np.zeros((specs[1]),dtype=np.int32)
	indptr=np.zeros((specs[0]+1),dtype=np.int32)

	# collect support
	_libShortCutSolver.Tools_Collect_VarList(
	        CClasses._toCInteger32Matrix(indices),CClasses._toCInteger32Matrix(indptr),
        	ct.c_long(pointerVarList[0]), ct.c_int(1))
	
	return (indices,indptr)
	


#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
# Modularized Setup Shielding Methods (Benchmark Variants)


def Setup_Shielding_Tree_Benchmark(NeighXindices,NeighXindptr,\
	lBottom, lTop, xpos, pointerY, pointerYpos, pointerYradii):

	pointer=np.zeros((1),dtype=np.int64)
	status=_libShortCutSolver.Setup_Shielding_Tree_Benchmark(\
		CClasses._toCInteger32Matrix(NeighXindices),CClasses._toCInteger32Matrix(NeighXindptr),\
		ct.c_int(lBottom), ct.c_int(lTop),
		CClasses._toCDoubleMatrix(xpos),
		ct.c_long(pointerY), CClasses._toCInteger64Matrix(pointerYpos), CClasses._toCInteger64Matrix(pointerYradii),
		CClasses._toCInteger64Matrix(pointer)\
		)

	result={"status" : status, "pointer" : pointer[0]}
	return result

def Setup_Shielding_TreePEucl_Benchmark(NeighXindices,NeighXindptr,\
	lBottom, lTop, xpos, pointerY, pointerYpos, pointerYradii, p, slack=0.):

	pointer=np.zeros((1),dtype=np.int64)
	status=_libShortCutSolver.Setup_Shielding_TreePEucl_Benchmark(\
		CClasses._toCInteger32Matrix(NeighXindices),CClasses._toCInteger32Matrix(NeighXindptr),\
		ct.c_int(lBottom), ct.c_int(lTop),\
		CClasses._toCDoubleMatrix(xpos),\
		ct.c_long(pointerY), CClasses._toCInteger64Matrix(pointerYpos), CClasses._toCInteger64Matrix(pointerYradii),\
		ct.c_double(p), ct.c_double(slack),\
		CClasses._toCInteger64Matrix(pointer)\
		)

	result={"status" : status, "pointer" : pointer[0]}
	return result

def Setup_Shielding_TreeSphere_Benchmark(NeighXindices, NeighXindptr,\
	lBottom, lTop, xpos, pointerY, pointerYpos, pointerYradii, p):

	pointer=np.zeros((1),dtype=np.int64)
	status=_libShortCutSolver.Setup_Shielding_TreeSphere_Benchmark(\
		CClasses._toCInteger32Matrix(NeighXindices),CClasses._toCInteger32Matrix(NeighXindptr),\
		ct.c_int(lBottom), ct.c_int(lTop),\
		CClasses._toCDoubleMatrix(xpos),\
		ct.c_long(pointerY), CClasses._toCInteger64Matrix(pointerYpos), CClasses._toCInteger64Matrix(pointerYradii),\
		ct.c_double(p),\
		CClasses._toCInteger64Matrix(pointer)\
		)

	result={"status" : status, "pointer" : pointer[0]}
	return result

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
# Interaction with active solver

def SolverStep(pointer,steps=1):
	status=_libShortCutSolver.ShortCutSolverStep(ct.c_long(pointer),steps)
	return {"status" : status}

def SolverGetReport(pointer):
	report=_libShortCutSolver.ShortCutSolverGetReport(ct.c_long(pointer))
	result={}
	for fname,ftype in report._fields_:
		result[fname]=eval("report."+fname)
	return result

def SolverClose(pointer):
	_libShortCutSolver.ShortCutSolverClose(ct.c_long(pointer))


def SolverGetSupport(pointer):
	# allocate small array for results of memory request call
	supportPointer=np.zeros((1),dtype=np.int64)
	supportSpecs=np.zeros((2),dtype=np.int32)
	# request how much memory is needed
	_libShortCutSolver.ShortCutSolverGetSupport_Request(ct.c_long(pointer),CClasses._toCInteger64Matrix(supportPointer),
		CClasses._toCInteger64Matrix(supportSpecs))
	
	# allocate space for varList
	varList=(\
		np.zeros((supportSpecs[1]),dtype=np.int32),
		np.zeros((supportSpecs[0]+1),dtype=np.int32)
		)
	
	# collect support
	_libShortCutSolver.Tools_Collect_VarList(
	        CClasses._toCInteger32Matrix(varList[0]),CClasses._toCInteger32Matrix(varList[1]),
        	ct.c_long(supportPointer[0]), ct.c_int(1))
	
	return varList

def SolverGetXVars(pointer):
	# allocate small array for results of memory request call
	varListPointer=np.zeros((1),dtype=np.int64)
	varListSpecs=np.zeros((2),dtype=np.int32)
	# request how much memory is needed
	_libShortCutSolver.ShortCutSolverGetXVars_Request(ct.c_long(pointer),CClasses._toCInteger64Matrix(varListPointer),
		CClasses._toCInteger64Matrix(varListSpecs))
	
	# allocate space for varList
	varList=(\
		np.zeros((varListSpecs[1]),dtype=np.int32),
		np.zeros((varListSpecs[0]+1),dtype=np.int32)
		)
	
	# collect support
	_libShortCutSolver.Tools_Collect_VarList(
	        CClasses._toCInteger32Matrix(varList[0]),CClasses._toCInteger32Matrix(varList[1]),
        	ct.c_long(varListPointer[0]), ct.c_int(0))
	
	return varList

########################################################################################################################
########################################################################################################################
########################################################################################################################
# Low level support methods


def _getFullVarList(xres, yres):
	indices=np.zeros((xres,yres),dtype=np.int32)
	indices[...]=np.arange(yres,dtype=np.int32)
	indices=indices.reshape((xres*yres))
	indPtr=np.arange(xres+1,dtype=np.int32)*yres
	return (indices,indPtr)


def GetGridNeighbours(dims,torusDim=0):
	neighbourPointer=np.zeros((1),dtype=np.int64)
	neighbourSpecs=np.zeros((2),dtype=np.int32)
	
	_libShortCutSolver.Tools_GetGridNeighbours_Request(\
		CClasses._toCInteger32Matrix(dims),CClasses._toCInteger64Matrix(neighbourPointer),\
		CClasses._toCInteger32Matrix(neighbourSpecs),\
		ct.c_int(torusDim))
		
	neighbourIndices=np.zeros((neighbourSpecs[1]),dtype=np.int32)
	neighbourIndptr=np.zeros((neighbourSpecs[0]+1),dtype=np.int32)

	_libShortCutSolver.Tools_Collect_VarList(
		CClasses._toCInteger32Matrix(neighbourIndices),CClasses._toCInteger32Matrix(neighbourIndptr),
		ct.c_long(neighbourPointer[0]), ct.c_int(1))
	
	return (neighbourIndices,neighbourIndptr)

########################################################################################################################
########################################################################################################################
# check shielding condition

def VerifyShielding(c,xVarsIndices,xVarsIndptr,neighboursIndices,neighboursIndptr,xMap):

	missesPointer=np.zeros((1),dtype=np.int64)
	missesSpecs=np.zeros((2),dtype=np.int32)

	result=_libShortCutSolver.Test_VerifyShielding(CClasses._toCDoubleMatrix(c),\
		CClasses._toCInteger32Matrix(xVarsIndices), CClasses._toCInteger32Matrix(xVarsIndptr),\
		CClasses._toCInteger32Matrix(neighboursIndices), CClasses._toCInteger32Matrix(neighboursIndptr),\
		CClasses._toCInteger32Matrix(xMap),\
		CClasses._toCInteger64Matrix(missesPointer), CClasses._toCInteger32Matrix(missesSpecs)\
		)

	missesIndices=np.zeros((missesSpecs[1]),dtype=np.int32)
	missesIndptr=np.zeros((missesSpecs[0]+1),dtype=np.int32)

	_libShortCutSolver.Tools_Collect_VarList(
		CClasses._toCInteger32Matrix(missesIndices),CClasses._toCInteger32Matrix(missesIndptr),
		ct.c_long(missesPointer[0]), ct.c_int(1))
	
	return (result,(missesIndices,missesIndptr))



def VerifyShieldingDuplex(c,xVarsIndices,xVarsIndptr,\
		neighboursXIndices,neighboursXIndptr,\
		neighboursYIndices,neighboursYIndptr,\
		xSupportIndices,xSupportIndptr):

	missesPointer=np.zeros((1),dtype=np.int64)
	missesSpecs=np.zeros((2),dtype=np.int32)

	result=_libShortCutSolver.Test_VerifyShieldingDuplex(CClasses._toCDoubleMatrix(c),\
		CClasses._toCInteger32Matrix(xVarsIndices), CClasses._toCInteger32Matrix(xVarsIndptr),\
		CClasses._toCInteger32Matrix(neighboursXIndices), CClasses._toCInteger32Matrix(neighboursXIndptr),\
		CClasses._toCInteger32Matrix(neighboursYIndices), CClasses._toCInteger32Matrix(neighboursYIndptr),\
		CClasses._toCInteger32Matrix(xSupportIndices), CClasses._toCInteger32Matrix(xSupportIndptr),\
		CClasses._toCInteger64Matrix(missesPointer), CClasses._toCInteger32Matrix(missesSpecs)\
		)
	
	missesIndices=np.zeros((missesSpecs[1]),dtype=np.int32)
	missesIndptr=np.zeros((missesSpecs[0]+1),dtype=np.int32)

	_libShortCutSolver.Tools_Collect_VarList(
		CClasses._toCInteger32Matrix(missesIndices),CClasses._toCInteger32Matrix(missesIndptr),
		ct.c_long(missesPointer[0]), ct.c_int(1))
	
	return (result,(missesIndices,missesIndptr))



def VerifyDualConstraints(c,alpha,beta,slack=1E-12):

	missesPointer=np.zeros((1),dtype=np.int64)
	missesSpecs=np.zeros((2),dtype=np.int32)


	result=_libShortCutSolver.Test_VerifyDualConstraints(\
		CClasses._toCDoubleMatrix(c),\
		CClasses._toCDoubleMatrix(alpha),\
		CClasses._toCDoubleMatrix(beta),\
		ct.c_double(slack),\
		CClasses._toCInteger64Matrix(missesPointer), CClasses._toCInteger32Matrix(missesSpecs)\
		)

	missesIndices=np.zeros((missesSpecs[1]),dtype=np.int32)
	missesIndptr=np.zeros((missesSpecs[0]+1),dtype=np.int32)

	_libShortCutSolver.Tools_Collect_VarList(
		CClasses._toCInteger32Matrix(missesIndices),CClasses._toCInteger32Matrix(missesIndptr),
		ct.c_long(missesPointer[0]), ct.c_int(1))
	
	return (result,(missesIndices,missesIndptr))

########################################################################################################################
########################################################################################################################
########################################################################################################################
# High level multi-scale solver

def MultiscaleSolver_SetupLevel(partitionX, partitionY, cList,\
	methodSetup_CouplingHandler, methodSetup_SubSolver, methodSetup_Shielding,\
	varList=None, nLayer=0,\
	Verbose=False, massKey1="mass", massKey2="mass", checkMethod=1):
	
	## prepare c,muX,muY,varList
	if Verbose:
		print("nLayer="+str(nLayer),flush=True)
		print("\tpreparing",flush=True)

	cGen=cList[nLayer]
	muXGen=partitionX.layers[nLayer][massKey1]
	muYGen=partitionY.layers[nLayer][massKey2]

	if varList is None:
		varListParam=_getFullVarList(muXGen.shape[0],muYGen.shape[0])
	else:
		varListParam=varList
	

	result={"status":0}

	## set up coupling handler
	if Verbose:
		print("\tsetting up coupling handler",flush=True)
	result_CouplingHandler=methodSetup_CouplingHandler(cGen,varListParam[0],varListParam[1])
	
	result["pointer_couplinghandler"]=result_CouplingHandler["pointer"]
	result["result_couplinghandler"]=result_CouplingHandler
	
	## set up sub-solver
	if Verbose:
		print("\tsetting up sub solver",flush=True)
	result_SubSolver=methodSetup_SubSolver(muXGen,muYGen,result["pointer_couplinghandler"])

	result["pointer_subsolver"]=result_SubSolver["pointer"]
	result["result_subsolver"]=result_SubSolver
	result["status"]=result_SubSolver["status"]
	if result["status"]!=0:
		return result

	## set up shield generator
	if Verbose:
		print("\tsetting up shielding generator",flush=True)
	result_Shielding=methodSetup_Shielding(nLayer,partitionX,partitionY)

	result["pointer_shielding"]=result_Shielding["pointer"]
	result["result_shielding"]=result_Shielding
	result["status"]=result_Shielding["status"]
	if result["status"]!=0:
		return result

	## set up short cut solver
	result_Solver=Setup_ShortCutSolver(result["pointer_couplinghandler"],\
		result["pointer_subsolver"],result["pointer_shielding"],checkMethod=checkMethod)

	result["pointer"]=result_Solver["pointer"]
	result["status"]=result_Solver["status"]
	if result["status"]!=0:
		return result

	
	return result	
	


def MultiscaleSolver(partitionX, partitionY, cList,\
	methodSetup_Refinement, methodSetup_CouplingHandler, methodSetup_SubSolver, methodSetup_Shielding,\
	nLayerFinal=None, nLayerInitial=0, varList=None,\
	Verbose=False,\
	maxSteps=100,  massKey1="mass", massKey2="mass", collectReports=False, measureTimes=False, checkMethod=1,\
	stepwiseAnalysis=False):
	
	# general initialization
	nLayer=nLayerInitial-1 # set to -1 such that first pseudo refinement gets it right
	if nLayerFinal is None:
		nLayerFinal=partitionX.nlayers-1


	if nLayerInitial>nLayerFinal:
		ValueError("Need nLayerInitial<=nLayerFinal")
		
	reportList=[]
	timeList=[]


	varListVar=varList;
	
	# iterating over finer scales
	while nLayer<nLayerFinal:
		nLayer=nLayer+1

		if measureTimes:
			time1=datetime.datetime.now()

		# do refinement (if necessary)
		if nLayer>nLayerInitial:
			# logging / benchmarking
			if Verbose:
				print("refining var list",flush=True)
			varListVar=methodSetup_Refinement(varListCoarse[0],varListCoarse[1],nLayer-1,nLayer-1)
		
		# initialization of solver
		result=MultiscaleSolver_SetupLevel(partitionX, partitionY, cList,\
			methodSetup_CouplingHandler, methodSetup_SubSolver, methodSetup_Shielding,\
			varList=varListVar,\
			nLayer=nLayer, Verbose=Verbose, massKey1=massKey1, massKey2=massKey2, checkMethod=checkMethod)
		pointer=result["pointer"]

		# logging
		if Verbose:
			print("\tsolving",flush=True)
					
		
		# do actual short-cut solving
		
		if stepwiseAnalysis:
			reportList.append([])
			steps=0
			status={"status":-1}
			while (steps<maxSteps) and (status["status"]!=0):
				status=SolverStep(pointer,1)
				report=SolverGetReport(pointer)
				reportList[-1].append(report)
				steps+=1
		else:
			status=SolverStep(pointer,maxSteps)
		
		
		# post-processing
		if status["status"]==0:
			if Verbose:
				print("\tsolved.",flush=True)
		else:
			print("Solving not successful. Final status: "+str(status))
			return result

		if measureTimes:
			time2=datetime.datetime.now()
			timeList.append((time2-time1).total_seconds())
	
		# get intermediate report
		if collectReports and (not stepwiseAnalysis):
			report=SolverGetReport(pointer)
			reportList.append(report)
			if Verbose:
				print(report,flush=True)
		
		if nLayer<nLayerFinal:
			# extract support (for refinement)
			varListCoarse=SolverGetSupport(pointer)
		
		if nLayer<nLayerFinal:
			if Verbose:
				print("close solver",flush=True)
			SolverClose(pointer)
	
	return (result,reportList,timeList)

########################################################################################################################
########################################################################################################################
########################################################################################################################
# Prepare setup methods

def getMethodSetup_Refinement(pointerX, pointerY, clib):
	return lambda vl0,vl1,lx,ly : clib.RefineVarList(pointerX,pointerY, vl0, vl1, lx, ly) 


# Coupling Handler

def getMethodSetup_CouplingHandler_SemiDense():
	return lambda c,VLXindices,VLXindptr :\
		Setup_CouplingHandler_SemiDense(c, VLXindices, VLXindptr, mu=None)

def getMethodSetup_CouplingHandler_Sparse_fullC():
	return lambda c,VLXindices,VLXindptr :\
		Setup_CouplingHandler_Sparse_fullC(c, VLXindices, VLXindptr)

# in this method we abuse the parameter c:
#	it will be a tuple containing (xPos,yPos)
def getMethodSetup_CouplingHandler_Sparse_dynamicC(method_setup_costFunctionProvider):
	def setup_coupling_handler(c,VLXindices, VLXindptr):
		CFPpointer=method_setup_costFunctionProvider(*c)
		xres=c[0].shape[0]
		yres=c[1].shape[0]
		result_CH=Setup_CouplingHandler_Sparse_dynamicC(VLXindices, VLXindptr,\
				xres, yres, CFPpointer)
		return {"status" :  result_CH["status"], "pointer" : result_CH["pointer"], "pointerCFP" : CFPpointer}
	return setup_coupling_handler


# Shielding

def getMethodSetup_Shielding_Grid():
	return lambda nLayer, partitionX, partitionY:\
		Setup_Shielding_Grid(partitionX.dims[nLayer], partitionY.dims[nLayer])

def getMethodSetup_Shielding_Padding(width=1):
	return lambda nLayer, partitionX, partitionY:\
		Setup_Shielding_Padding(partitionX.dims[nLayer], partitionY.dims[nLayer],width)


def getMethodSetup_Shielding_Tree(neighboursList,pointerY,pointerYpos,pointerYradii):
	return lambda nLayer, partitionX, partitionY:\
		Setup_Shielding_Tree(neighboursList[nLayer][0],neighboursList[nLayer][1],\
			nLayer, 0, partitionX.layers[nLayer]["pos"],\
			pointerY, pointerYpos, pointerYradii)

def getMethodSetup_Shielding_Tree_Torus(neighboursList,pointerY,pointerYpos,pointerYradii,pointerYTorusRadii,torusRadii,torusDim):
	return lambda nLayer, partitionX, partitionY:\
		Setup_Shielding_Tree_Torus(neighboursList[nLayer][0],neighboursList[nLayer][1],\
			nLayer, 0, partitionX.layers[nLayer]["pos"],\
			pointerY, pointerYpos, pointerYradii,pointerYTorusRadii,torusRadii,torusDim)

def getMethodSetup_Shielding_TreeNoise(neighboursList,pointerY,pointerYpos,pointerYradii,pointerCList, c_eta_list, c_lambda_list):
	return lambda nLayer, partitionX, partitionY:\
		Setup_Shielding_TreeNoise(neighboursList[nLayer][0],neighboursList[nLayer][1],\
			nLayer, 0, partitionX.layers[nLayer]["pos"],\
			pointerY, pointerYpos, pointerYradii,\
			pointerCList, c_eta_list[nLayer], c_lambda_list[nLayer]\
			)

def getMethodSetup_Shielding_TreePEucl(neighboursList,pointerY,pointerYpos,pointerYradii,p,slack=0.):
	return lambda nLayer, partitionX, partitionY:\
		Setup_Shielding_TreePEucl(neighboursList[nLayer][0],neighboursList[nLayer][1],\
			nLayer, 0, partitionX.layers[nLayer]["pos"],\
			pointerY, pointerYpos, pointerYradii,p,slack=slack)

def getMethodSetup_Shielding_TreeSphere(neighboursList,pointerY,pointerYpos,pointerYradii,p):
	return lambda nLayer, partitionX, partitionY:\
		Setup_Shielding_TreeSphere(neighboursList[nLayer][0],neighboursList[nLayer][1],\
			nLayer, 0, partitionX.layers[nLayer]["pos"],\
			pointerY, pointerYpos, pointerYradii,p)

def getMethodSetup_Shielding_TreeReflector(neighboursList,pointerY,pointerYpos,pointerYradii):
	return lambda nLayer, partitionX, partitionY:\
		Setup_Shielding_TreeReflector(neighboursList[nLayer][0],neighboursList[nLayer][1],\
			nLayer, 0, partitionX.layers[nLayer]["pos"],\
			pointerY, pointerYpos, pointerYradii)


# Shielding Benchmark

def getMethodSetup_Shielding_Tree_Benchmark(neighboursList,pointerY,pointerYpos,pointerYradii):
	return lambda nLayer, partitionX, partitionY:\
		Setup_Shielding_Tree_Benchmark(neighboursList[nLayer][0],neighboursList[nLayer][1],\
			nLayer, 0, partitionX.layers[nLayer]["pos"],\
			pointerY, pointerYpos, pointerYradii)

def getMethodSetup_Shielding_TreePEucl_Benchmark(neighboursList,pointerY,pointerYpos,pointerYradii,p,slack=0.):
	return lambda nLayer, partitionX, partitionY:\
		Setup_Shielding_TreePEucl_Benchmark(neighboursList[nLayer][0],neighboursList[nLayer][1],\
			nLayer, 0, partitionX.layers[nLayer]["pos"],\
			pointerY, pointerYpos, pointerYradii,p,slack=slack)

def getMethodSetup_Shielding_TreeSphere_Benchmark(neighboursList,pointerY,pointerYpos,pointerYradii,p):
	return lambda nLayer, partitionX, partitionY:\
		Setup_Shielding_TreeSphere_Benchmark(neighboursList[nLayer][0],neighboursList[nLayer][1],\
			nLayer, 0, partitionX.layers[nLayer]["pos"],\
			pointerY, pointerYpos, pointerYradii,p)
