import numpy as np
from ..Common import *
import warnings

def getChildrenPosList(parentPos):
	relpos=getPoslistNCube([2 for i in range(len(parentPos))],dtype=np.int32)
	offset=np.array([2*p for p in parentPos])
	result=[offset+pos for pos in relpos]
	return result

def getChildrenLists(pointList, indexList, level, parentPos, box):
	childrenPosList=getChildrenPosList(parentPos)
	availableList=np.full(indexList.shape,1,dtype=np.bool)
	childrenLists=[]
	posLists=[]
	for pos in childrenPosList:
		inbox=isInBox(box[0]+(box[1]-box[0])*pos/(2**(level+1)),box[0]+(box[1]-box[0])*(pos+1)/(2**(level+1)), pointList)
		inbox=np.logical_and(inbox,availableList)
		availableList=np.logical_and(availableList,np.invert(inbox))
		preChildrenList=np.nonzero(inbox)[0]
		
		childrenLists.append(indexList[preChildrenList])
		posLists.append(pos)
		
	return (childrenLists,posLists)


def getDims(levels,dim,final):
	"""Get grid dimension arrays for hierarchical partition on a grid. For the upper layers it's simply
	2**(nLevels-1-level) along each axis, lowest level is given by given full fine grid."""
	result=np.zeros((levels,dim),dtype=np.int32)
	for i in range(0,levels-1):
		result[i]=2**i
	result[-1]=final
	return result
	
	
def getSignalPointer(partition, signame, lBottom=None):
	if lBottom is None:
		lBottom=partition.nlayers-1

	result=np.zeros((lBottom+1),dtype=np.int64)
	for i in range(lBottom+1):
		result[i]=partition.layers[i][signame].ctypes.data
	return result

def getListPointer(signal):
	result=np.zeros((len(signal)),dtype=np.int64)
	for i in range(len(signal)):
		result[i]=signal[i].ctypes.data
	return result

# Two constants to describe two different ways to enumerate children in THierarchicalPartition
# Tree is adaptive: child nodes will only be created if elements are in there, and enumerated such that children of one node have
#	successive numbers.
# Grid is more rigid and grid based: All child nodes will be created, regardless of whether they are empty. They will be enumrated according
#	to their position on the grid.

# The former is in general more efficient and dynamically adapts to point clouds in higher dimensions, which model lower dimensional
#	sub-manifolds.
# The latter is needed for simplified shielding neighbourhood construction on R^n

THPMode_Tree=0
THPMode_Grid=1

class THierarchicalPartition:
	def __init__(self,_points,\
		box=None,boxThresh=1E-10,childMode=THPMode_Tree):
		
		self.points=_points
		if box is not None:
			self._box=box
		else:
			# automatically build max bounding box.
			self._box=[np.amin(self.points,axis=0)-boxThresh,np.amax(self.points,axis=0)+boxThresh]
		
		self._childMode=childMode
		
		self.npoints=self.points.shape[0]
		self.ndim=self.points.shape[1]
		
	def reset(self):
		self.layers=[{"posCode" : np.zeros((1,self.ndim),dtype=np.int32),\
			"children" : [np.array([],dtype=np.int32)], "leaves" : [np.arange(self.npoints,dtype=np.int32)],\
			"parent" : np.array([-1],dtype=np.int32)}]
		self.nlayers=1
		self.cardLayers=[1]

	def setup(self, nMiddleLayers, finestDims=None, finestDimsWarning=True):
		self.reset()
		for i in range(nMiddleLayers):
			self.refine()
		self.addAtomicLayer()
		
		# if in THPMode_Grid, add grid dimensions at each level
		if self._childMode==THPMode_Grid:
			if finestDims is None:
				if finestDimsWarning:
					warnings.warn("No finestDims specified for HierarchicalPartition in HTPMode_Grid. \
						Cannot run ShortCut solver.")
			else:
				self.dims=getDims(self.nlayers,self.ndim,finestDims)
				
		    
		
	def refineCell(self,layer,cellid):
		"""For a given layer cell (specified by layer & cellid), returns two lists.
		childrenLists and posLists. First giving leaves for each child, second gives the poscode."""
		indices=self.layers[layer]["leaves"][cellid]
		pos=self.layers[layer]["posCode"][cellid]
		childrenLists=getChildrenLists(self.points[indices],indices,layer,pos,self._box)
		return childrenLists

	def refine(self):
		"""Add new refined layer at bottom of Hierarchical Partition."""
		
		if self._childMode==THPMode_Tree:
			# create new empty dummy layer
			newLayer={"posCode" : [], \
				"children" : [], "leaves" : [], "parent" : []}
			# offset is needed to assign each new child a unique id
			offset=0

		elif self._childMode==THPMode_Grid:
			nChildNodes=2**(self.nlayers*self.ndim)
			# create grid-based dummy layer
			newLayer={"posCode" : np.zeros((nChildNodes,self.ndim),dtype=np.int32),\
				"children" : [[] for i in range(nChildNodes)],\
				"leaves" : [[] for i in range(nChildNodes)],\
				"parent": np.zeros((nChildNodes),dtype=np.int32)\
				}
				
			# needed to transform posCode into location in childList
			posVector=np.array([2**(self.nlayers*(self.ndim-i-1)) for i in range(self.ndim)],dtype=np.int32)
		
		# go through each cell in latest layer
		layer=self.nlayers-1
		for cellid in range(self.cardLayers[layer]):
			# get refinement data of element
			childrenLists,posLists=self.refineCell(layer,cellid)
			
			if self._childMode==THPMode_Tree:
				# count proper children (len(cl)>0):
				childrenCount=0
				# go through new child elements and add them to new layer
				for cl,pl in zip(childrenLists,posLists):
					if len(cl)>0:
						childrenCount+=1
						newLayer["posCode"].append(pl)
						#newLayer["mass"].append(np.sum(self.mu[cl]))
						newLayer["children"].append(np.array([],dtype=np.int32))
						newLayer["leaves"].append(cl)
						newLayer["parent"].append(cellid)
				# set children list in parent on last layer
				self.layers[layer]["children"][cellid]=np.arange(childrenCount,dtype=np.int32)+offset
				# increase offset by children in last cell
				offset+=childrenCount
			
			elif self._childMode==THPMode_Grid:
				# reset children list of future parent (may already be np.array, but we need dynamical appending)
				self.layers[layer]["children"][cellid]=[]
				# go through new child elements and add them to new layer
				for cl,pl in zip(childrenLists,posLists):
					if len(cl)==0:
						warnings.warn("Node without children occured.")
					childPos=np.sum(posVector*pl)
					newLayer["posCode"][childPos]=pl
					newLayer["leaves"][childPos]=cl
					newLayer["parent"][childPos]=cellid
					# add to children list in parent on last layer
					self.layers[layer]["children"][cellid].append(childPos)
				
				# transform children list back to np.array
				self.layers[layer]["children"][cellid]=np.array(self.layers[layer]["children"][cellid],dtype=np.int32)
				
		
		if self._childMode==THPMode_Tree:
			# transform some of the data to np.array
			for ind,dtyp in zip(["posCode","parent"],[np.int32,np.int32]):
				newLayer[ind]=np.array(newLayer[ind],dtype=dtyp)
		
		# add layer to list, update nlayers, cardLayers
		self.layers.append(newLayer)
		self.nlayers+=1
		self.cardLayers.append(newLayer["parent"].shape[0])
	    
	def addAtomicLayer(self):
		"""Add a final layer to Hierarchical Partition in which each element is an atomic singleton of one leave.
		This layer does not provide the full specifications of upper layers.
		poscode, children and leaves do not make sense here anymore."""
		# create layer object. pos, mass are trivial on this layer.
		newLayer={\
			"parent" : np.full((self.points.shape[0]),-1,dtype=np.int32)}
		# parent needs to be set by looking at last layer
		layer=self.nlayers-1
		for cellid in range(self.cardLayers[layer]):
			# at this stage, the children are equal to leaves
			self.layers[layer]["children"][cellid]=self.layers[layer]["leaves"][cellid]
			# set parent value of all children
			newLayer["parent"][self.layers[layer]["children"][cellid]]=cellid
		
		# add layer to list, update nlayers, cardLayers
		self.layers.append(newLayer)
		self.nlayers+=1
		self.cardLayers.append(self.points.shape[0])
	

	def getSignalPos(self):
		signal=[None for i in range(self.nlayers)]
		signal[self.nlayers-1]=self.points.astype(dtype=np.double)
		for i in range(0,self.nlayers-1):
			signal[i]=self._box[0]+(self._box[1]-self._box[0])*(self.layers[i]["posCode"]+0.5)/(2**(i))
		
		return signal
	
	def getSignalRadii(self,axis=None):
		if axis is None:
			_axis=slice(0,None)
		else:
			_axis=axis
			
		# compute diameter of whole box (of selected axis)
		diagonal=self._box[1][_axis]-self._box[0][_axis]
		diagonal=np.sqrt(np.sum(diagonal*diagonal))
		signal=[None for i in range(self.nlayers)]
		for i in range(0,self.nlayers-1):
			signal[i]=np.full((self.cardLayers[i]),diagonal/(2**(i+1)),dtype=np.double)
		
		return signal
		
		
	def addSignal(self,signame,signal):
		for i in range(self.nlayers):
			self.layers[i][signame]=signal[i]
			
	def addSignalPos(self):
		self.addSignal("pos",self.getSignalPos())
		
	def addSignalRadii(self):
		self.addSignal("radii",self.getSignalRadii())

##################################################################################################################

def GetPartition(pos, depth, childMode=THPMode_Tree, finestDims=None, box=None,\
	signal_pos=False, signal_radii=False, clib=None, export=False, mu=None, massKey="mass",\
	verbose=False, finestDimsWarning=True):
	"""Convenience routine to do many standard things for setting up a hierarchical partition."""
    
	partition=THierarchicalPartition(pos,childMode=childMode,box=box)
	partition.setup(depth, finestDims=finestDims, finestDimsWarning=finestDimsWarning)
    
	if verbose:
		print(partition.cardLayers)
    
	if signal_pos:
		partition.addSignalPos()
    
	if signal_radii:
		partition.addSignalRadii()
    
	if (clib is not None) and export:
		pointer=clib.Export(partition)
    
	if (mu is not None) and (clib is not None) and export:
		partition.addSignal(massKey,clib.GetSignalMass(pointer,partition,mu))

	if (clib is not None) and export:
		return (partition,pointer)
    
	return partition

def GetHierarchicalCost(partitionX, partitionY, func):
	return [func(partitionX.layers[i]["pos"],partitionY.layers[i]["pos"]) for i in range(partitionX.nlayers)]

def GetNearestNeighbours(metric, nNeighbours):
	neighbours=np.argsort(metric,axis=1)[:,1:1+nNeighbours]
	neighboursIndsep=np.arange(0,neighbours.shape[0]+1,dtype=np.int32)*neighbours.shape[1]
	neighboursIndices=neighbours.astype(dtype=np.int32).flatten()
	return (neighboursIndices,neighboursIndsep)

def GetHierarchicalNearestNeighbours(metricList, nNeighbours,layerTop=1):
	neighboursList=[None for i in range(layerTop)]
	for i in range(layerTop,len(metricList)):
		neighboursList.append(GetNearestNeighbours(metricList[i],nNeighbours))
	return neighboursList
	
	
###############################################################################################
# Pseudo layers for gradual refinement

def AddPseudoLayer(partition,nLayer,addPseudoRadii=False):
	oldLayer=partition.layers[nLayer]
	# create basic coarse and fine layer objects
	newLayerCoarse={\
			"children" : [np.array([i],dtype=np.int32) for i in range(partition.cardLayers[nLayer])],\
			"parent" : oldLayer["parent"]}
	newLayerFine={\
			"parent" : np.arange(partition.cardLayers[nLayer],dtype=np.int32)}
    
	if "children" in oldLayer.keys():
		newLayerFine["children"]=oldLayer["children"]
    
	if "leaves" in oldLayer.keys():
		newLayerCoarse["leaves"]=oldLayer["leaves"]
		newLayerFine["leaves"]=oldLayer["leaves"]
	else:
		# need to invent pseudo leaves layer for the new coarse layer, just copy children layer
		newLayerCoarse["leaves"]=newLayerCoarse["children"]
        
    
	# add additional signals/objects if required
	for signal in ["pos","radii","posCode"]:
		if signal in oldLayer.keys():
			newLayerCoarse[signal]=oldLayer[signal]
			newLayerFine[signal]=oldLayer[signal]
        
	# add pseudo radii for new coarse layer if required
	if addPseudoRadii and (oldLayer["radii"] is None):
		newLayerCoarse["radii"]=np.zeros((partition.cardLayers[nLayer]),dtype=np.double)
    
	# restructure partition object
	partition.layers = partition.layers[0:nLayer] + [newLayerCoarse,newLayerFine] + partition.layers[nLayer+1:]
	partition.nlayers+=1
	partition.cardLayers = partition.cardLayers[0:nLayer] +\
			[partition.cardLayers[nLayer],partition.cardLayers[nLayer]]+\
			partition.cardLayers[nLayer+1:]
    
def AddGradualRefinement(partitionList,nLayer,addPseudoRadii=False):
	"""Adds a suitable set of pseudo-layers to the partitions in partitionList to create a gradual refinement
	from nLayer to nLayer+1."""
	for i in range(1,len(partitionList)):
		for j in range(len(partitionList)):
			if j<i:
				AddPseudoLayer(partitionList[j],nLayer+i,addPseudoRadii=addPseudoRadii)
			else:
				AddPseudoLayer(partitionList[j],nLayer,addPseudoRadii=addPseudoRadii)
