import numpy as np
import scipy
from . import Common
import OTTools

def Iterate(kernelRows,kernelCols,u,v,mu,nu,n=1):
	for i in range(n):
		v[...]=nu/kernelCols.dot(u)
		u[...]=mu/kernelRows.dot(v)


def GetMarginals(kernelRows,kernelCols,u,v):
	m1=u*(kernelRows.dot(v))
	m2=v*(kernelCols.dot(u))
	return [m1,m2]
        

def ErrorMarginLInf(kernelRows,kernelCols,u,v,mu,nu):
	marginals=GetMarginals(kernelRows,kernelCols,u,v)
	result=0
	for m,n in zip(marginals,[mu,nu]):
		result=max(result,np.max(np.abs(m-n)))
	return result

def ErrorMarginL1(kernelRows,kernelCols,u,v,mu,nu):
	marginals=GetMarginals(kernelRows,kernelCols,u,v)
	result=0
	for m,n in zip(marginals,[mu,nu]):
		result+=np.sum(np.abs(m-n))
	return result


def ThreshParticles(particles,minMass):
	keepPos=np.nonzero((particles[2]>minMass))
	particles=[p[keepPos] for p in particles]
	return particles

def GetParticles(pi,minMass=None):
	posArray=OTTools.GetFullPosArray(pi.data,pi.indices,pi.indptr)
	particles=(posArray[:,0],posArray[:,1],pi.data)
	
	if minMass is not None:
		particles=ThreshParticles(particles,minMass)
	
	return particles


def interpolateEuclidean(particles,xpos,ypos,t):
	"""particles: (x0_id,x1_id,m0)"""
	return ((1-t)*xpos[particles[0]]+t*ypos[particles[1]],particles[2])

def interpolateDistances(particles,t):
	"""Gives Wasserstein interpolation of a list of Dirac particles.
	
	Arguments:
	particles= lists of (dist,m)
	
	Return list of (x,m):
	x: fraction of distances travelled,
	m: intermediate masses.
	
	For standard Wasserstein space this function is trivial, it exists for interface consistence with other models."""

	
	return (np.full(particles[1].shape,t,dtype=np.double),particles[1])

