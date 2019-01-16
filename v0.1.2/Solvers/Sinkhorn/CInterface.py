import os
import numpy as np
import scipy
import ctypes as ct
import CClasses

import OTTools.HierarchicalPartition as HierarchicalPartition

_package_directory=os.path.dirname(os.path.abspath(__file__))

_libSinkhorn = ct.cdll.LoadLibrary(_package_directory+'/Release/libSinkhorn.so')

# C++ multivar interface

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
# Basic Aux Functions

def CollectPosVarList(pointer, specs=None, doDelete=True):
	"""Collects a TKernelHandler_Pos from c++. pointer gives address in memory, specs is an array with [total,dim].
		doDelete specifies whether to free c++ object upon retrieval."""
	
	# delete parameter
	if doDelete:
		deletePar=1
	else:
		deletePar=0

	if specs is None:
		specsParam=GetSpecsPosVarList(pointer)
	else:
		specsParam=specs

	# allocate space for varList
	varList=(\
		np.zeros((specsParam[0]),dtype=np.double),\
		np.zeros((specsParam[0],specsParam[1]),dtype=np.int32)\
		)\
	
	status=_libSinkhorn.Tools_Collect_PosVarList(
			CClasses._toCDoubleMatrix(varList[0]),
			CClasses._toCInteger32Matrix(varList[1]),
			ct.c_long(pointer), ct.c_int(deletePar))
	
	if status!=0:
		raise ValueError("An error occurred during retrieval of varList: "+str(status))
	
	return varList


def ExportPosVarList(varList):
	pointer=np.zeros((1),dtype=np.int64)
	status= _libSinkhorn.Tools_Export_PosVarList(\
		CClasses._toCDoubleMatrix(varList[0]), CClasses._toCInteger32Matrix(varList[1]),\
		CClasses._toCInteger64Matrix(pointer),
		)
	if status!=0:
		raise ValueError("An error occurred during exporting of varList: "+str(status))
	return pointer[0]

def DeletePosVarList(pointer):
	return _libSinkhorn.Tools_Delete_PosVarList(\
		ct.c_long(pointer)\
		)

def GetSpecsPosVarList(pointer):
	specs=np.zeros((2),dtype=np.int32)
	_libSinkhorn.Tools_GetSpecs_PosVarList(\
		ct.c_long(pointer),\
		CClasses._toCInteger32Matrix(specs)\
		)
	return specs

###################################


def GetMarginalPos(pointerKernel, pointerListScaling, res, axis):
	mu=np.zeros((res[axis]),dtype=np.double)
	_libSinkhorn.Tools_GetMarginal(ct.c_long(pointerKernel), CClasses._toCInteger64Matrix(pointerListScaling),\
		CClasses._toCDoubleMatrix(mu), CClasses._toCInteger32Matrix(res), ct.c_uint(axis))
	return mu


def GetMarginalsPos(pointerKernel,pointerListScaling,res):
	
	return [GetMarginalPos(pointerKernel, pointerListScaling, res, a) \
			for a in range(len(pointerListScaling))]



def Kernel_Scale(pointerKernel, pointerListScaling):
	_libSinkhorn.Kernel_Scale(ct.c_long(pointerKernel), CClasses._toCInteger64Matrix(pointerListScaling))


def Kernel_Exponentiate(pointerKernel, eps):
	_libSinkhorn.Kernel_Exponentiate(ct.c_long(pointerKernel), ct.c_double(eps))

##############################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
# CostFunctionProviderSetup

def _getPseudoDualList(alpha, lTop):
	"""Aux function if one wants to set up a cost function provider which is only evaluated at finest scale and thus only
		finest alpha is required. So this creates an empty pointer list, only setting the finest level properly."""
	result=np.zeros((lTop+1),dtype=np.int64)
	result[lTop]=alpha.ctypes.data
	return result

#######################
# SquaredEuclidean

#def Setup_CostFunctionProvider_SquaredEuclidean(pointerPosX,pointerPosY,posDim, lBottom,\
#		pointerRadiiX=None,pointerRadiiY=None,pointerAlpha=None,pointerBeta=None,\
#		alphaFinest=None, betaFinest=None):
#	
#	pointerPos=np.zeros((2),dtype=np.int64)
#	pointerPos[0]=pointerPosX.ctypes.data
#	pointerPos[1]=pointerPosY.ctypes.data
#	
#	pointerRadii=np.zeros((2),dtype=np.int64)
#	if pointerRadiiX is not None:
#		pointerRadii[0]=pointerRadiiX.ctypes.data
#		pointerRadii[1]=pointerRadiiY.ctypes.data
#	else:
#		pointerRadii[0]=0
#		pointerRadii[0]=0
#	
#	pointerDuals=np.zeros((2),dtype=np.int64)
#	pointerPseudoDuals=None
#	if pointerAlpha is not None:
#		pointerDuals[0]=pointerAlpha.ctypes.data
#		pointerDuals[1]=pointerBeta.ctypes.data
#	elif alphaFinest is not None:
#		# only need pseudo lists where only entry lBottom is valid. other levels will not be called.
#		pointerPseudoDuals=[_getPseudoDualList(fine,lBottom) for fine in [alphaFinest,betaFinest]]
#		pointerDuals[0]=pointerPseudoDuals[0].ctypes.data
#		pointerDuals[1]=pointerPseudoDuals[1].ctypes.data
#	else:
#		pointerDuals[0]=0
#		pointerDuals[1]=0
#	
#	# return CFP type value as well, for convenience, also return all allocated arrays to stop automatic garbage collection
#	return ((_libSinkhorn.Setup_CostFunctionProvider_SquaredEuclidean(
#		CClasses._toCInteger64Matrix(pointerPos),CClasses._toCInteger64Matrix(pointerRadii),CClasses._toCInteger64Matrix(pointerDuals),
#		ct.c_int(posDim), ct.c_int(lBottom)),1),\
#		(pointerPos,pointerRadii,pointerDuals,pointerPseudoDuals))

#######################
# SquaredEuclidean

def Setup_CostFunctionProvider_SquaredEuclidean(pointerPos,posDim, lBottom,\
		pointerRadii=None,pointerAlpha=None,\
		alphaFinest=None, weight=1.):
	
	if pointerRadii is not None:
		pointerRadiiParam=pointerRadii
	else:
		pointerRadiiParam=np.zeros_like(pointerPos)

	pointerPseudoDuals=None
	if pointerAlpha is not None:
		pointerAlphaParam=pointerAlpha
	elif alphaFinest is not None:
		# only need pseudo lists where only entry lBottom is valid. other levels will not be called.
		pointerPseudoDuals=[_getPseudoDualList(fine,lBottom) for fine in alphaFinest]
		pointerAlphaParam=np.array([ptr.ctypes.data for ptr in pointerPseudoDuals],dtype=np.int64)
	else:
		pointerAlphaParam=np.zeros_like(pointerPos)
	
	
	# return CFP type value as well, for convenience. also return all allocated arrays to stop automatic garbage collection
	return ((
			_libSinkhorn.Setup_CostFunctionProvider_SquaredEuclidean(
					CClasses._toCInteger64Matrix(pointerPos),CClasses._toCInteger64Matrix(pointerRadiiParam),\
					CClasses._toCInteger64Matrix(pointerAlphaParam),
					ct.c_int(posDim), ct.c_int(lBottom), ct.c_double(weight)),
			1),\
			(pointerRadiiParam,pointerAlphaParam,pointerPseudoDuals))


#######################
# ColorTransport: Squared Euclidean + RGB

def Setup_CostFunctionProvider_Color_SquaredEuclidean_RGB(pointerPos,posDim, lBottom, lTop, colorWeight,\
		pointerRadii=None,pointerAlpha=None,\
		alphaFinest=None,\
		FR_mode=False, FR_kappa=1.):
	
	if pointerRadii is not None:
		pointerRadiiParam=pointerRadii
	else:
		pointerRadiiParam=np.zeros_like(pointerPos)

	pointerPseudoDuals=None
	if pointerAlpha is not None:
		pointerAlphaParam=pointerAlpha
	elif alphaFinest is not None:
		# only need pseudo lists where only entry lBottom is valid. other levels will not be called.
		pointerPseudoDuals=[_getPseudoDualList(fine,lBottom) for fine in alphaFinest]
		pointerAlphaParam=np.array([ptr.ctypes.data for ptr in pointerPseudoDuals],dtype=np.int64)
	else:
		pointerAlphaParam=np.zeros_like(pointerPos)
		
	
	# FR mode
	if FR_mode:
		_FR_mode=1
	else:
		_FR_mode=0
	
	
	# return CFP type value as well, for convenience, also return all allocated arrays to stop automatic garbage collection
	return ((
			_libSinkhorn.Setup_CostFunctionProvider_Color_SquaredEuclidean_RGB(
				CClasses._toCInteger64Matrix(pointerPos),CClasses._toCInteger64Matrix(pointerRadiiParam),\
				CClasses._toCInteger64Matrix(pointerAlphaParam),
				ct.c_int(posDim), ct.c_int(lBottom),\
				ct.c_int(lTop),ct.c_double(colorWeight),ct.c_int(_FR_mode),ct.c_double(FR_kappa)),\
			4),\
			(pointerRadiiParam,pointerAlphaParam,pointerPseudoDuals)
			)


#######################
# ColorTransport: Squared Euclidean + HSV_HS

def Setup_CostFunctionProvider_Color_SquaredEuclidean_HSV_HS(pointerPos,posDim, lBottom, liftMode, colorWeight,\
		pointerRadiiPos=None, pointerRadiiHue=None, pointerRadiiVal=None,\
		pointerAlpha=None,\
		alphaFinest=None,\
		FR_mode=False, FR_kappa=1.):
	
	if pointerRadiiPos is not None:
		pointerRadiiPosParam=pointerRadiiPos
	else:
		pointerRadiiPosParam=np.zeros_like(pointerPos)

	if pointerRadiiHue is not None:
		pointerRadiiHueParam=pointerRadiiHue
	else:
		pointerRadiiHueParam=np.zeros_like(pointerPos)

	if pointerRadiiVal is not None:
		pointerRadiiValParam=pointerRadiiVal
	else:
		pointerRadiiValParam=np.zeros_like(pointerPos)


	pointerPseudoDuals=None
	if pointerAlpha is not None:
		pointerAlphaParam=pointerAlpha
	elif alphaFinest is not None:
		# only need pseudo lists where only entry lBottom is valid. other levels will not be called.
		pointerPseudoDuals=[_getPseudoDualList(fine,lBottom) for fine in alphaFinest]
		pointerAlphaParam=np.array([ptr.ctypes.data for ptr in pointerPseudoDuals],dtype=np.int64)
	else:
		pointerAlphaParam=np.zeros_like(pointerPos)
		
	
	# FR mode
	if FR_mode:
		_FR_mode=1
	else:
		_FR_mode=0
	
	
	# return CFP type value as well, for convenience, also return all allocated arrays to stop automatic garbage collection
	return ((
			_libSinkhorn.Setup_CostFunctionProvider_Color_SquaredEuclidean_HSV_HS(
				CClasses._toCInteger64Matrix(pointerPos),\
				CClasses._toCInteger64Matrix(pointerRadiiPosParam),\
				CClasses._toCInteger64Matrix(pointerRadiiHueParam),\
				CClasses._toCInteger64Matrix(pointerRadiiValParam),\
				CClasses._toCInteger64Matrix(pointerAlphaParam),
				ct.c_int(posDim), ct.c_int(lBottom),\
				ct.c_int(liftMode),ct.c_double(colorWeight),ct.c_int(_FR_mode),ct.c_double(FR_kappa)),\
			4),\
			(pointerRadiiPosParam,pointerRadiiHueParam,pointerRadiiValParam,pointerAlphaParam,pointerPseudoDuals)
			)


#######################
# SquaredEuclidean + FR

def Setup_CostFunctionProvider_SquaredEuclideanWF(pointerPos,posDim, lBottom,\
		pointerRadii=None,pointerAlpha=None,\
		alphaFinest=None, FR_kappa=1., FR_cMax=1E10):
	
	if pointerRadii is not None:
		pointerRadiiParam=pointerRadii
	else:
		pointerRadiiParam=np.zeros_like(pointerPos)

	pointerPseudoDuals=None
	if pointerAlpha is not None:
		pointerAlphaParam=pointerAlpha
	elif alphaFinest is not None:
		# only need pseudo lists where only entry lBottom is valid. other levels will not be called.
		pointerPseudoDuals=[_getPseudoDualList(fine,lBottom) for fine in alphaFinest]
		pointerAlphaParam=np.array([ptr.ctypes.data for ptr in pointerPseudoDuals],dtype=np.int64)
	else:
		pointerAlphaParam=np.zeros_like(pointerPos)
	
	
	# return CFP type value as well, for convenience, also return all allocated arrays to stop automatic garbage collection
	return ((\
			_libSinkhorn.Setup_CostFunctionProvider_SquaredEuclideanWF(
					CClasses._toCInteger64Matrix(pointerPos),CClasses._toCInteger64Matrix(pointerRadiiParam),\
					CClasses._toCInteger64Matrix(pointerAlphaParam),
					ct.c_int(posDim), ct.c_int(lBottom), ct.c_double(FR_kappa), ct.c_double(FR_cMax)\
			),\
			3),\
			(pointerRadiiParam,pointerAlphaParam,pointerPseudoDuals)\
			)

#######################
# SquaredEuclideanBarycenter

def Setup_CostFunctionProvider_SquaredEuclideanBarycenter(pointerPos,weightList, posDim, lBottom,pointerRadii=None,pointerAlpha=None,
	alphaFinest=None):
	
	if pointerRadii is not None:
		pointerRadiiParam=pointerRadii
	else:
		pointerRadiiParam=np.zeros_like(pointerPos)

	pointerPseudoDuals=None
	if pointerAlpha is not None:
		pointerAlphaParam=pointerAlpha
	elif alphaFinest is not None:
		# only need pseudo lists where only entry lBottom is valid. other levels will not be called.
		pointerPseudoDuals=[_getPseudoDualList(fine,lBottom) for fine in alphaFinest]
		pointerAlphaParam=np.array([ptr.ctypes.data for ptr in pointerPseudoDuals],dtype=np.int64)
	else:
		pointerAlphaParam=np.zeros_like(pointerPos)
	
	# return CFP type value as well, for convenience
	return ((_libSinkhorn.Setup_CostFunctionProvider_SquaredEuclideanBarycenter(
		CClasses._toCInteger64Matrix(pointerPos),CClasses._toCInteger64Matrix(pointerRadiiParam),\
		CClasses._toCInteger64Matrix(pointerAlphaParam),CClasses._toCDoubleMatrix(weightList),\
		ct.c_int(posDim), ct.c_int(lBottom)),2),(pointerRadiiParam,pointerAlphaParam,pointerPseudoDuals))


#######################
# Coulomb Cost

def Setup_CostFunctionProvider_Coulomb(pointerPos, charges, posDim, lBottom,pointerRadii=None,pointerAlpha=None,
	alphaFinest=None):
	
	if pointerRadii is not None:
		pointerRadiiParam=pointerRadii
	else:
		pointerRadiiParam=np.zeros_like(pointerPos)

	pointerPseudoDuals=None
	if pointerAlpha is not None:
		pointerAlphaParam=pointerAlpha
	elif alphaFinest is not None:
		# only need pseudo lists where only entry lBottom is valid. other levels will not be called.
		pointerPseudoDuals=[_getPseudoDualList(fine,lBottom) for fine in alphaFinest]
		pointerAlphaParam=np.array([ptr.ctypes.data for ptr in pointerPseudoDuals],dtype=np.int64)
	else:
		pointerAlphaParam=np.zeros_like(pointerPos)
	
	# return CFP type value as well, for convenience
	return ((
			_libSinkhorn.Setup_CostFunctionProvider_Coulomb(
					CClasses._toCInteger64Matrix(pointerPos),CClasses._toCInteger64Matrix(pointerRadiiParam),\
					CClasses._toCInteger64Matrix(pointerAlphaParam),\
					ct.c_int(posDim), ct.c_int(lBottom),CClasses._toCDoubleMatrix(charges)),\
			0),\
		(pointerRadiiParam,pointerAlphaParam,pointerPseudoDuals))


#######################
# Spherical Reflector

def Setup_CostFunctionProvider_Reflector_Spherical(pointerPos,posDim, lBottom,\
		pointerRadii=None,pointerAlpha=None,\
		alphaFinest=None):
	
	if pointerRadii is not None:
		pointerRadiiParam=pointerRadii
	else:
		pointerRadiiParam=np.zeros_like(pointerPos)

	pointerPseudoDuals=None
	if pointerAlpha is not None:
		pointerAlphaParam=pointerAlpha
	elif alphaFinest is not None:
		# only need pseudo lists where only entry lBottom is valid. other levels will not be called.
		pointerPseudoDuals=[_getPseudoDualList(fine,lBottom) for fine in alphaFinest]
		pointerAlphaParam=np.array([ptr.ctypes.data for ptr in pointerPseudoDuals],dtype=np.int64)
	else:
		pointerAlphaParam=np.zeros_like(pointerPos)
	
	
	# return CFP type value as well, for convenience. also return all allocated arrays to stop automatic garbage collection
	return ((
			_libSinkhorn.Setup_CostFunctionProvider_Reflector_Spherical(
					CClasses._toCInteger64Matrix(pointerPos),CClasses._toCInteger64Matrix(pointerRadiiParam),\
					CClasses._toCInteger64Matrix(pointerAlphaParam),
					ct.c_int(posDim), ct.c_int(lBottom)),
			1),\
			(pointerRadiiParam,pointerAlphaParam,pointerPseudoDuals))


#######################
# Interpolator


def Setup_CostFunctionProvider_Interpolator(coarse, fine, pointerListPartition, q, lBottom ,pointerAlpha=None,
	alphaFinest=None):

	pointerPseudoDuals=None
	if pointerAlpha is not None:
		pointerAlphaParam=pointerAlpha
	elif alphaFinest is not None:
		# only need pseudo lists where only entry lBottom is valid. other levels will not be called.
		pointerPseudoDuals=[_getPseudoDualList(fine,lBottom) for fine in alphaFinest]
		pointerAlphaParam=np.array([ptr.ctypes.data for ptr in pointerPseudoDuals],dtype=np.int64)
	else:
		pointerAlphaParam=np.zeros_like(pointerListPartition)
	

	# return CFP type value as well, for convenience
	return ((_libSinkhorn.Setup_CostFunctionProvider_Interpolator(\
		ct.c_long(coarse[0][0]), ct.c_long(fine[0][0]),
		CClasses._toCInteger64Matrix(pointerListPartition),
		ct.c_double(q),
		CClasses._toCInteger64Matrix(pointerAlphaParam)),\
		3),(coarse,fine,pointerAlphaParam,pointerPseudoDuals))

###########################################################################################
# Some Aux Functions for Interaction with CostFunctionProviders

def Delete_CostFunctionProvider(pointerCFP):
	_libSinkhorn.Tools_Delete_CostFunctionProvider(ct.c_long(pointerCFP))


def Evaluate_CostFunctionProvider(CFPAddr, pos, layer):
	c=np.zeros((pos.shape[0]),dtype=np.double)
	status=_libSinkhorn.Tools_Evaluate_CostFunctionProvider(
		ct.c_long(CFPAddr), CClasses._toCInteger32Matrix(pos), CClasses._toCDoubleMatrix(c), ct.c_int(layer))
	if status is not 0:
		raise ValueError("Error during evaluation of CostFunctionProvider.")
	return c


#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
# High Level Functions

def CheckDualConstraints_Pos_generateKernel(pointerPartition, pointerCFP,\
	slack, lBottom=None):
	
	if lBottom is None:
		lBottomParam=len(alphaList)-1
	else:
		lBottomParam=lBottom

	nProbs=len(pointerPartition)
	
	specs=np.zeros((2),dtype=np.int32)
	varListPointer=np.zeros((1),dtype=np.int64)
	
	
	_libSinkhorn.Check_DualConstraints_Pos(CClasses._toCInteger64Matrix(pointerPartition),
		ct.c_int(lBottomParam),
		ct.c_double(slack),
		ct.c_long(pointerCFP),
		CClasses._toCInteger64Matrix(specs), CClasses._toCInteger64Matrix(varListPointer))
	
	return (varListPointer,specs)


def CheckDualConstraints_Pos(pointerPartition, pointerCFP,\
	slack, lBottom=None):
	
	(pointer,specs)=CheckDualConstraints_Pos_generateKernel(pointerPartition, pointerCFP,\
		slack, lBottom)
	return CollectPosVarList(pointer,specs,doDelete=True)


def RefineVarList_CSR_Pos_generateKernel(pointerPartition, pointerCFP, lTop, varList):
	
	specs=np.zeros((2),dtype=np.int32)
	varListPointer=np.zeros((1),dtype=np.int64)
	

	_libSinkhorn.Refine_VarList_CSR_Pos(CClasses._toCInteger64Matrix(pointerPartition),
		ct.c_int(lTop),
		ct.c_long(pointerCFP),
		CClasses._toCInteger32Matrix(varList[0]), CClasses._toCInteger32Matrix(varList[1]),
		CClasses._toCInteger64Matrix(specs), CClasses._toCInteger64Matrix(varListPointer))
		
	return (varListPointer,specs)


def RefineVarList_CSR_Pos(pointerPartition, pointerCFP, lTop, varList):

	
	(pointer,specs)=RefineVarList_CSR_Pos_generateKernel(pointerPartition, pointerCFP, lTop, varList)
	return CollectPosVarList(pointer,specs,doDelete=True)


def RefineVarList_Pos_Pos_generateKernel(pointerPartition, pointerCFP, lTop, pointerKernel):
	
	specs=np.zeros((2),dtype=np.int32)
	varListPointer=np.zeros((1),dtype=np.int64)
	

	_libSinkhorn.Refine_VarList_Pos_Pos(CClasses._toCInteger64Matrix(pointerPartition),\
		ct.c_int(lTop),\
		ct.c_long(pointerCFP),\
		ct.c_long(pointerKernel),\
		CClasses._toCInteger64Matrix(specs), CClasses._toCInteger64Matrix(varListPointer))
		
	return (varListPointer,specs)



def ReEvaluateVarList_generateVarLists(pointerPartition, pointerCFP, lTop, varList):
	
	specs=np.zeros((2),dtype=np.int32)
	varListPointer=np.zeros((1),dtype=np.int64)
	

	_libSinkhorn.ReEvaluate_VarList_Pos(\
		CClasses._toCInteger64Matrix(pointerPartition),\
		ct.c_int(lTop),\
		ct.c_long(pointerCFP),\
		CClasses._toCInteger32Matrix(varList[0]), CClasses._toCInteger32Matrix(varList[1]),\
		CClasses._toCInteger64Matrix(specs), CClasses._toCInteger64Matrix(varListPointer))
	
	return (varListPointer,specs)


def ReEvaluateVarList(pointerPartition, pointerCFP, lTop, varList):
	
	(pointer,specs)=ReEvaluateVarList_generateVarLists(pointerPartition, pointerCFP, lTop, varList)
	return CollectPosVarList(pointer,specs,doDelete=True)




def GetDenseCosts_Pos_generateKernel(pointerPartition, pointerCFP,\
	lBottom=None):
	
	if lBottom is None:
		lBottomParam=len(alphaList)-1
	else:
		lBottomParam=lBottom
	
	specs=np.zeros((2),dtype=np.int32)
	varListPointer=np.zeros((1),dtype=np.int64)
			
	_libSinkhorn.Tools_GetDenseCosts_Pos(CClasses._toCInteger64Matrix(pointerPartition),
		ct.c_int(lBottomParam),
		ct.c_long(pointerCFP),
		CClasses._toCInteger64Matrix(specs), CClasses._toCInteger64Matrix(varListPointer))
	
	return (varListPointer,specs)


def GetDenseCosts_Pos(pointerPartition, pointerCFP,\
	lBottom=None):
	
	(pointer,specs)=GetDenseCosts_Pos_generateKernel(pointerPartition, pointerCFP,\
		lBottom)
	return CollectPosVarList(pointer,specs,doDelete=True)


def Iterate(pointerKernel, pointerListScaling, pointerListMu, res, n):
	_libSinkhorn.Iterate(\
		ct.c_long(pointerKernel),\
		CClasses._toCInteger64Matrix(pointerListScaling),\
		CClasses._toCInteger64Matrix(pointerListMu),\
		CClasses._toCInteger32Matrix(res),\
		ct.c_int(n))


