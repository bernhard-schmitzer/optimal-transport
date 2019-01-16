import numpy as np
import scipy
from . import Common
from .. import Aux
import OTTools
from . import OT

def interpolateEuclidean(particles,xpos,ypos,t,M0,M1,zeroThresh=1E-10):
	"""particles: (x0_id,x1_id,m0)"""
	x0=xpos[particles[0]]
	x1=ypos[particles[1]]
	m0=particles[2]
	deltaX=np.linalg.norm(x0-x1,axis=1)

	(x,m)=interpolateDistances((deltaX,m0),t,M0,M1,zeroThresh=zeroThresh)
	x=x0+np.einsum(x,[0],x1-x0,[0,1],[0,1])
	return (x,m)


def interpolateDistances(particles,t,M0,M1,zeroThresh=1E-10):
	"""Gives limit FR_kappa->infty Wasserstein-Fisher-Rao interpolation of a list of Dirac particles.
	
	Arguments:
	particles= lists of (dist,m0), between probability normalized measures
	M0, M1: total masses of original marginal measures
	
	Return list of (x,m):
	x: fraction of distances travelled,
	m: intermediate masses."""

	# if masses are numerically equal, just call standard OT interpolation
	if np.abs(M0-M1)<=zeroThresh:
		(x,m)=OT.interpolateDistances(particles,t)
		m=M0*m
		return (x,m)

	deltaX=particles[0]	
	m0=particles[1]
	
	if t==0:
		return (np.full(m0.shape,0.,dtype=np.double),M0*m0)
	if t==1:
		return (np.full(m0.shape,1.,dtype=np.double),M1*m0)
        
	# compute aux constants
	A=M0+M1-2*np.sqrt(M0*M1)
	B=M0-np.sqrt(M0*M1)

	# interpolate mass
	m=((1-t)*np.sqrt(M0)+t*np.sqrt(M1))**2*m0


	# interpolate location
	# location factor
	x=np.sqrt(M0*M1)*A*t/((B-A*t)*B)
	# constant array
	x=np.full(m0.shape,x,dtype=np.double)

    	# return list of intermediate positions and masses at time t
	return (x,m)
