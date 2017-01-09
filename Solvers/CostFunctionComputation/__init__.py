import os
import numpy as np
import ctypes as ct
import CClasses
import OTTools.HierarchicalPartition.CInterface as OTHPCInterface
import OTTools.HierarchicalPartition as HierarchicalPartition

_package_directory=os.path.dirname(os.path.abspath(__file__))


_libCFC = ct.cdll.LoadLibrary(_package_directory+'/Release/libCostFunctionComputation.so')


def Export(partition):
	"""Export a hierarchical partition to c++. Return pointer."""
	return OTHPCInterface.Export(_libCFC, partition)


def Close(pointer):
	"""Clean up an exported hiearchical partition."""
	return OTHPCInterface.Close(_libCFC, pointer)


def GetSignalMass(pointer, partition, mu):
	"""Compute hierarchical multiscale representation of measure mu."""
	return OTHPCInterface.GetSignalMass(_libCFC, pointer, partition, mu)

	
def RefineVarList(pointerX, pointerY, VLXindices, VLXindptr, nLayerX, nLayerY):
	"""For a sparse varlist (VLXindices,VLXindptr) on layers (nLayerX,nLayerY), do refinement by one layer."""
	# allocate small array for results of memory request call
	requiredVariables=np.zeros((2),dtype=np.int32)

	# request how much memory is needed
	_libCFC.refineVarListRequestMemory(ct.c_long(pointerX),ct.c_long(pointerY),\
		CClasses._toCInteger32Matrix(VLXindices),CClasses._toCInteger32Matrix(VLXindptr),\
		ct.c_int(nLayerX),ct.c_int(nLayerY),\
		CClasses._toCInteger32Matrix(requiredVariables))
	
	# allocate space for varList
	varListFine=(\
		np.zeros((requiredVariables[1]),dtype=np.int32),
		np.zeros((requiredVariables[0]+1),dtype=np.int32)
		)

	# request refined var list	
	_libCFC.refineVarList(ct.c_long(pointerX),ct.c_long(pointerY),\
		CClasses._toCInteger32Matrix(VLXindices),CClasses._toCInteger32Matrix(VLXindptr),\
		ct.c_int(nLayerX),ct.c_int(nLayerY),\
		CClasses._toCInteger32Matrix(varListFine[0]),CClasses._toCInteger32Matrix(varListFine[1])\
		)
	
	return varListFine


def RefineSignal(partitionX, pointerX, f, lTop):
	"""Compute refinement of signal f living on layer lTop."""

	
	# allocate space for fine signal
	signal=np.zeros((partitionX.cardLayers[lTop+1]),dtype=np.double)
	
	# call clib routine	
	_libCFC.refine_signal(ct.c_long(pointerX),
		CClasses._toCDoubleMatrix(f), CClasses._toCDoubleMatrix(signal), ct.c_int(lTop))
		
	return signal

MODE_MAX=1
MODE_MIN=0


def GetSignalMinList(partitionX, pointerX, f, mode=MODE_MIN, lTop=None, lBottom=None):
	"""Compute hierarchical representation of signal f on layer lBottom to all coarser layers. Modes: MODE_MIN, MODE_MAX."""

	if lTop is None:
		lTopParam=0
	else:
		lTopParam=lTop
	
	if lBottom is None:
		lBottomParam=partitionX.nlayers-1
	else:
		lBottomParam=lBottom

	
	# allocate space for hierarchical signal
	signal=[np.zeros((partitionX.cardLayers[layer]),dtype=np.double) for layer in range(lTopParam,lBottomParam)]
	# add finest layer at end
	signal=signal+[f]
	
	# set up list of pointers to layers
	signalPointer=np.zeros((lBottomParam-lTopParam+1),dtype=np.int64)
	for i in range(lBottomParam-lTopParam+1):
		signalPointer[i]=signal[i].ctypes.data
	
	# call clib routine	
	_libCFC.propagate_signalFunction(ct.c_long(pointerX),
		CClasses._toCInteger64Matrix(signalPointer), ct.c_int(lTopParam), ct.c_int(lBottomParam), ct.c_int(mode))
		
	return signal

def GetCMinList(partitionX, partitionY, pointerX, pointerY, c, mode=MODE_MIN, lTop=None, lBottom=None):
	"""Compute hierarchical representation of cost signal living over product space of two hierarchical partitions. Modes: MODE_MIN, MODE_MAX."""
	if lTop is None:
		lTopParam=0
	else:
		lTopParam=lTop
	
	if lBottom is None:
		lBottomParam=partitionX.nlayers-1
	else:
		lBottomParam=lBottom

	
	# allocate space for hierarchical signal
	signal=[np.zeros((partitionX.cardLayers[layer],partitionY.cardLayers[layer]),dtype=np.double) for layer in range(lTopParam,lBottomParam)]
	# add finest layer at end
	signal=signal+[c]

	# set up list of pointers to layers
	signalPointer=np.zeros((lBottomParam-lTopParam+1),dtype=np.int64)
	for i in range(lBottomParam-lTopParam+1):
		signalPointer[i]=signal[i].ctypes.data

			
	_libCFC.propagate_costFunction(ct.c_long(pointerX),ct.c_long(pointerY),
		CClasses._toCInteger64Matrix(signalPointer), ct.c_int(lTopParam), ct.c_int(lBottomParam), ct.c_int(mode))
		
	return signal


def CheckDualConstraints(pointerX, pointerY, cList, alphaList, betaList, slack=None, slackX=None, slackY=None, lTop=None, lBottom=None):
	if lTop is None:
		lTopParam=0
	else:
		lTopParam=lTop
	
	if lBottom is None:
		lBottomParam=len(cList)-1
	else:
		lBottomParam=lBottom

	# determine whether to run adaptive or standard version
	if slack is not None:
		slackMode=0
	elif (slackX is not None) and (slackY is not None):
		slackMode=1
		slackXPointer=HierarchicalPartition.getListPointer(slackX)
		slackYPointer=HierarchicalPartition.getListPointer(slackY)
	else:
		raise ValueError("Either slack or (slackX,slackY) must be provided.")

	cPointer=HierarchicalPartition.getListPointer(cList)
	alphaPointer=HierarchicalPartition.getListPointer(alphaList)
	betaPointer=HierarchicalPartition.getListPointer(betaList)
	
	specs=np.zeros((2),dtype=np.int32)
	varListPointer=np.zeros((1),dtype=np.int64)
	
	if slackMode==0:
		_libCFC.check_dualConstraints(ct.c_long(pointerX),ct.c_long(pointerY),
			CClasses._toCInteger64Matrix(cPointer),
			CClasses._toCInteger64Matrix(alphaPointer), CClasses._toCInteger64Matrix(betaPointer),
			ct.c_int(lTopParam), ct.c_int(lBottomParam), ct.c_double(slack),
			CClasses._toCInteger64Matrix(specs), CClasses._toCInteger64Matrix(varListPointer))
	elif slackMode==1:
		_libCFC.check_dualConstraints_adaptive(ct.c_long(pointerX),ct.c_long(pointerY),
			CClasses._toCInteger64Matrix(cPointer),
			CClasses._toCInteger64Matrix(alphaPointer), CClasses._toCInteger64Matrix(betaPointer),
			ct.c_int(lTopParam), ct.c_int(lBottomParam),
			CClasses._toCInteger64Matrix(slackXPointer), CClasses._toCInteger64Matrix(slackYPointer),
			CClasses._toCInteger64Matrix(specs), CClasses._toCInteger64Matrix(varListPointer))
	
	# allocate space for varList
	varList=(\
		np.zeros((specs[1]),dtype=np.int32),
		np.zeros((specs[0]+1),dtype=np.int32)
		)
	
	# collect support
	_libCFC.Tools_Collect_VarList(
	        CClasses._toCInteger32Matrix(varList[0]),CClasses._toCInteger32Matrix(varList[1]),
        	ct.c_long(varListPointer[0]), ct.c_int(1))

	return varList

def GetCEffectiveValues(c, alpha, beta, indices, indptr):
	"""For entries specified by (indices,indptr) (in sparse.csr format) extract c[x,y]-alpha[x]-beta[y] as data array."""
	
	data=np.zeros(indices.shape,dtype=np.double)
		
	_libCFC.get_CEffectiveValues(
		CClasses._toCInteger32Matrix(indices), CClasses._toCInteger32Matrix(indptr),
		CClasses._toCDoubleMatrix(data),
		CClasses._toCDoubleMatrix(c), CClasses._toCDoubleMatrix(alpha), CClasses._toCDoubleMatrix(beta)\
		)
	
	return data
