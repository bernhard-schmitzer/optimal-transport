import numpy as np

import ctypes as ct
import CClasses

def Export(lib, partition):
	"""Exports a THierarchicalPartition object to a loaded c++ library via the interface specified in THierarchicalPartitionInterface.h"""
	pointerArray=np.zeros((1),dtype=np.int64)
	
	# Create new HierarchicalPartition object in the library
	res=lib.HierarchicalPartitionsCreate(ct.c_int(partition.nlayers),ct.c_int(partition.ndim),CClasses._toCInteger64Matrix(pointerArray))
	pointer=pointerArray[0];
	
	# Now go through list of layers and export them
	for layerNr in range(partition.nlayers):
		res=lib.HierarchicalPartitionsCreateLayer(ct.c_long(pointer), ct.c_int(layerNr),\
			ct.c_int(partition.cardLayers[layerNr]),\
			CClasses._toCInteger32Matrix(partition.layers[layerNr]["parent"])\
			);
		
		
		# if not in finest = atomic layer
		if layerNr<partition.nlayers-1:
			# now go through list of cells and export all children and leaves lists
			for cellNr in range(partition.cardLayers[layerNr]):
				lib.HierarchicalPartitionsCreateLayerAddSublists(ct.c_long(pointer), ct.c_int(layerNr),\
					ct.c_int(cellNr),\
					CClasses._toCInteger32Matrix(partition.layers[layerNr]["children"][cellNr]),\
					CClasses._toCInteger32Matrix(partition.layers[layerNr]["leaves"][cellNr]))

	return pointer

def GetSignalMass(lib, pointer, partition, mu):
	"""For a measure mu compute its hierarchical representation. Requires that partition is already exported to c++."""
	# allocate space for hierarchical mass array
	signal=[np.zeros((card),dtype=np.double) for card in partition.cardLayers]
	# set up list of pointers to layers
	signalPointer=np.zeros((partition.nlayers),dtype=np.int64)
	for i in range(partition.nlayers):
		signalPointer[i]=signal[i].ctypes.data
	
	# call clib routine	
	lib.HierarchicalPartitionsGetSignalMass(ct.c_long(pointer),CClasses._toCDoubleMatrix(mu),CClasses._toCInteger64Matrix(signalPointer))
	return signal
	

def Close(lib, pointer):
	"""Clean up exported partition."""
	lib.HierarchicalPartitionsClose(ct.c_long(pointer))

