import numpy as np
from .Common import *

def processDensity_Grid(x,totalMass=None,constOffset=None,keepZero=True,zeroThresh=1E-14):

	# process actual density
	
	# copy, cast to double and reshape
	img=x.astype(dtype=np.double).copy()
	shape=img.shape
	nPoints=np.prod(shape)
	dim=len(shape)
	img=img.reshape((nPoints))
	
	processDensity(img,totalMass=totalMass,constOffset=constOffset)

	# get grid pos information
	posList=getPoslistNCube(shape,dtype=np.int)
	posList=posList.reshape((nPoints,dim))

	# if desired, throw away points with zero mass
	if not keepZero:
		nonZeroPos=np.nonzero(img>zeroThresh)
		img=img[nonZeroPos]
		posList=posList[nonZeroPos]
		
		# if necessary, rescale mass once more
		processDensity(img, totalMass=totalMass, constOffset=None)

	return (img,posList)

def processDensity(x, totalMass=None, constOffset=None):
	# re-normalize and add offset if required
	if totalMass is not None:
		x[:]=totalMass*x/np.sum(x)
		if constOffset is not None:
			x[:]=x+constOffset
			x[:]=totalMass*x/np.sum(x)
	else:
		if constOffset is not None:
			x[:]=x+constOffset
			#x[:]=x/np.sum(x)

def TruncateMeasures(muX,muY,measures_scale):
	muX[...]=np.trunc(muX/measures_scale)
	muY[...]=np.trunc(muY/measures_scale)
	muXsum=muX.sum()
	muYsum=muY.sum()
	if muXsum<muYsum:
		muX[-1]+=muYsum-muXsum
	else:
		muY[-1]+=muXsum-muYsum

	if np.min(muX)<=0.:
		raise ValueError("muX min is <= 0.")
	if np.min(muY)<=0.:
		raise ValueError("muY min is <= 0.")
	if np.sum(muX)!=np.sum(muY):
		raise ValueError("Unequal masses.")


def getEuclideanCostFunction(x1,x2,p=2.):
	c=-2*np.einsum(x1,[0,2],x2,[1,2],[0,1])
	x1Sqr=np.sum(x1*x1,axis=1)
	x2Sqr=np.sum(x2*x2,axis=1)
	for i in range(x1.shape[0]):
		c[i,:]+=x1Sqr[i]+x2Sqr
	return np.power(c,p/2.)
	
def getSparseFromMask(mask,c=None):
	"""Extracts sparse index data from the array mask. It is assumed that mask has shape (xres,yres) and is boolean,
	True indicating that a given entry is selected and False otherwise.
	The result is (data,indices,indptr) as used by scipy.sparse.csr_matrix.
	If c is None, only (indices,indptr) are returned."""
	
	xres=mask.shape[0]
	total=np.sum(mask)
	# if c is given, reserve space for data
	if c is not None:
		extractData=True
		data=np.zeros((total),dtype=np.double)
	else:
		extractData=False
	
	# allocate space for indices and indptr
	indices=np.zeros((total),dtype=np.int32)
	indptr=np.zeros((xres+1),dtype=np.int32)
	
	# go through rows of mask / cost function
	offset=0
	for x in range(xres):
		currentIndices=np.nonzero(mask[x])[0]
		rowLen=currentIndices.shape[0]
		indices[offset:offset+rowLen]=currentIndices
		if extractData:
			data[offset:offset+rowLen]=c[x,currentIndices]
		
		offset+=rowLen
		indptr[x+1]=offset
	if extractData:
		return (data,indices,indptr)
	else:
		return (indices,indptr)

def getFullVarList(xres, yres):
	indices=np.zeros((xres,yres),dtype=np.int32)
	indices[...]=np.arange(yres,dtype=np.int32)
	indices=indices.reshape((xres*yres))
	indPtr=np.arange(xres+1,dtype=np.int32)*yres
	return (indices,indPtr)


