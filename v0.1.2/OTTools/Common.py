import numpy as np
import scipy

def getPoslistNCube(shape,dtype=np.double):
	"""Create list of positions in an n-dimensional cube of size shape."""
	ndim=len(shape)

	axGrids=[np.arange(i,dtype=dtype) for i in shape]
	prePos=np.array(np.meshgrid(*axGrids,indexing='ij'),dtype=dtype )
	# the first dimension of prepos is the dimension of the  posvector, the successive dimensions are in the cube
	# so need to move first axis to end, and then flatten
	pos=np.rollaxis(prePos,0,ndim+1)
	# flattening
	newshape=(-1,ndim)
	return (pos.reshape(newshape)).copy()

def GetProductGrid(g1,g2):
	"""For two point lists g1, g2 with shapes (M,m), (N,n) (e.g. M points in m dimensions), compute the product grid of shape (M*N,m+n)."""
	M,m=g1.shape
	N,n=g2.shape

	result=np.zeros((M,N,m+n),dtype=np.double)
	result[:,:,m:]=g2

	resultT=result.transpose((1,0,2))
	resultT[:,:,:m]=g1

	result=result.reshape((M*N,m+n))
	return result


def isInBox(x0,x1,y):
	"""Checks if the points in y are in the hypercube spanned by corners x0 and x1. (x0<=y<x1 to be precise.)
	Returns list of bool."""
	return np.logical_and( np.all(x0<=y,axis=1), np.all(y<=x1,axis=1) )
	
def GetFullPosArray(data,indices,indptr):
	"""In a standard sparse kernel rep, the indices only contain coords of axis [1,..,dim-1].
	This method returns a new indices array where the first axis is added. Useful for post-processing."""

	if indices.ndim==1:
		posDim=2
	else:
		posDim=indices.shape[1]+1
    
	newIndices=np.zeros((indices.shape[0],posDim),dtype=np.int32)

	newIndices[:,1:]=indices.reshape((indices.shape[0],posDim-1))
	for i in range(indptr.shape[0]-1):
		newIndices[indptr[i]:indptr[i+1],0]=i
    
	return newIndices
	
def GetBarycenterPoints(posArray,weightList,posList):
	posDim=posList[0].shape[1]
	result=np.zeros((posArray.shape[0],posDim),dtype=np.double)
	for i in range(len(weightList)):
		result[...]+=weightList[i]*posList[i][posArray[:,i]]
	return result
	
def GetGaussian1D(n,x0,sigma,totalMass=1.):
	"""Creates a 1d Gaussian, image size n, centered at x0, std dev sigma, normalized to totalMass. If None, density will be 1 at x0."""
	img=np.arange(n)
	img=np.exp(-0.5*((img-x0)/(sigma))**2)
	if totalMass is not None:
		img[...]=totalMass*img/np.sum(img)
	return img

def GetGaussian2D(n,x0,sigma,totalMass=1.):
	"""Creates a 2d Gaussian, image size n=(n[0],n[1]), centered at x0, std dev (sigma[0],sigma[1]) along first and second axis,
		normalized to totalMass. If None, density will be 1 at x0."""
	img0=GetGaussian1D(n[0],x0[0],sigma[0],totalMass=totalMass)
	img1=GetGaussian1D(n[1],x0[1],sigma[1],totalMass=totalMass)
	img=np.einsum(img0,[0],img1,[1],[0,1])
    
	if totalMass is not None:
		img[...]=totalMass*img/np.sum(img)
	return img


def ProjectInterpolation1D(pos,masses,xres):
	"""Projects list of particles with positions pos and masses masses onto 1D pixel grid of length xres,
	by spreading mass linearly between two closest pixels."""	
	indices=np.zeros((pos.shape[0],2),dtype=np.double)
	data=np.zeros((pos.shape[0],2),dtype=np.double)
	# set x indices
	indices[:,0]=np.floor(pos[:,0])
	indices[:,1]=np.floor(pos[:,0])+1
    
	# set data
	for i in range(2):
		data[:,i]=masses*(1-np.abs(pos[:,0]-indices[:,i]))
    
	# truncate index range
	indices[:]=np.maximum(indices,0)
	indices[:]=np.minimum(indices,xres-1)
    
	result=scipy.sparse.coo_matrix((data.flatten(),(indices.flatten(),np.zeros((pos.shape[0]*2),dtype=np.int))),\
			shape=(xres,1))
	return result

def ProjectInterpolation2D(pos,masses,xres,yres):
	"""Projects list of particles with positions pos and masses masses onto 2D pixel grid of shape (xres,yres),
	by spreading mass linearly between four closest pixels."""	
	indices=np.zeros((2,pos.shape[0],4),dtype=np.double)
	data=np.zeros((pos.shape[0],4),dtype=np.double)
	# set x indices
	indices[0,:,0]=np.floor(pos[:,0])
	indices[0,:,1]=np.floor(pos[:,0])+1
	indices[0,:,2]=np.floor(pos[:,0])
	indices[0,:,3]=np.floor(pos[:,0])+1
	# set y indices
	indices[1,:,0]=np.floor(pos[:,1])
	indices[1,:,1]=np.floor(pos[:,1])
	indices[1,:,2]=np.floor(pos[:,1])+1
	indices[1,:,3]=np.floor(pos[:,1])+1
	# set data
	for i in range(4):
		data[:,i]=masses*(1-np.abs(pos[:,0]-indices[0,:,i]))*(1-np.abs(pos[:,1]-indices[1,:,i]))

	# truncate index range
	indices[0,...]=np.maximum(indices[0,...],0)
	indices[0,...]=np.minimum(indices[0,...],xres-1)
	indices[1,...]=np.maximum(indices[1,...],0)
	indices[1,...]=np.minimum(indices[1,...],yres-1)

        
	result=scipy.sparse.coo_matrix((data.flatten(),(indices[0].flatten(),indices[1].flatten())),\
			shape=(xres,yres))
	return result

def AverageMongeMap(muX,pi,posY):
	"""Compute the center of mass of each row of coupling pi (scipy.sparse.csr_matrix), first marginal given by muX,
		positions of column points given by posY."""
	dim=posY.shape[1]
	xres=muX.shape[0]
	t=np.zeros((xres,dim),dtype=np.double)
	for i in range(dim):
		t[:,i]=pi.dot(posY[:,i])/muX
	return t
