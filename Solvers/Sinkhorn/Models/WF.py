import numpy as np
import scipy
from . import Common
from .. import Aux
import OTTools

def Iterate(kernelRows,kernelCols,alpha,beta,u,v,mu,nu,eps,FR_kappa,n):
	prefac=2*FR_kappa**2
	for k in range(n):
		conv=kernelCols.dot(u)
		if np.min(conv)==0:
			# maybe think of a way to act on zero cols as well.
			#posZero=np.nonzero(conv==0)
			posNonZero=np.nonzero(conv)
			v[posNonZero]=(nu[posNonZero]/conv[posNonZero])**(prefac/(prefac+eps))*np.exp(-beta[posNonZero]/(prefac+eps))
		else:
			v[:]=(nu/conv)**(prefac/(prefac+eps))*np.exp(-beta/(prefac+eps))
        
		conv=kernelRows.dot(v)
		if np.min(conv)==0:
			posNonZero=np.nonzero(conv)
			u[posNonZero]=(mu[posNonZero]/conv[posNonZero])**(prefac/(prefac+eps))*np.exp(-alpha[posNonZero]/(prefac+eps))
		else:
			u[:]=(mu/conv)**(prefac/(prefac+eps))*np.exp(-alpha/(prefac+eps))



def ScorePrimal(kernelRows,kernelCols,alpha,beta,u,v,mu,nu,eps,FR_kappa):

	prefac=2*FR_kappa**2
	
	alphaEff=alpha+eps*np.log(u)
	betaEff=beta+eps*np.log(v)
	m0=u*(kernelRows.dot(v))
	m1=v*(kernelCols.dot(u))

	return np.sum(m0*alphaEff)+np.sum(m1*betaEff)+prefac*(Common.ScoreKL(m0,mu)+Common.ScoreKL(m1,nu))


def ScoreDual(kernelRows,kernelCols,alpha,beta,u,v,mu,nu,eps,FR_kappa):

	prefac=2*FR_kappa**2

	alphaEff=alpha+eps*np.log(u)
	betaEff=beta+eps*np.log(v)

	convX=kernelRows.dot(v)
	convY=kernelCols.dot(u)

	result=0

	# colums with zero conv should effectively have dual variables +inf, i.e. they will make no contribution
	# in the following terms
	posNonZero=np.nonzero(convX)
	result+=prefac*(np.sum(mu)-np.sum(np.exp(-alphaEff[posNonZero]/prefac)*mu[posNonZero]))

	posNonZero=np.nonzero(convY)
	result+=prefac*(np.sum(nu)-np.sum(np.exp(-betaEff[posNonZero]/prefac)*nu[posNonZero]))


	return result

def ScorePDGap(kernelRows,kernelCols,alpha,beta,u,v,mu,nu,eps,FR_kappa):
	return ScorePrimal(kernelRows,kernelCols,alpha,beta,u,v,mu,nu,eps,FR_kappa)\
		-ScoreDual(kernelRows,kernelCols,alpha,beta,u,v,mu,nu,eps,FR_kappa)


def GetSemiCouplingScaling(mu,m0):
    s0=np.ones(mu.shape,dtype=np.double)
    nonZeroPos=np.nonzero(m0)
    s0[nonZeroPos]=mu[nonZeroPos]/m0[nonZeroPos]
    return s0


def GetSemiCouplings(pi,mu,nu,m0=None,m1=None):
	"""Get two semi-couplings from a soft-marginal coupling. Useful for post-processing.
		If marginals of pi are already known, can be given as m0, m1. If not given, they will
		be computed.
		Returns (pi0,pi1,sMin)
		Both semi-couplings are divided by minimal scaling factor, sMin, which is also returned.
		So actual semi couplings are given by pi0*sMin, pi1*sMin.
		This keeps the sparse module from restructuring the matrices and thus both keep same entries,
		which simplifies creating particle list."""
		
	if (m0 is None) or (m1 is None):
		_m0,_m1=Common.GetMarginals(pi)
	else:
		_m0,_m1=m0,m1

	# scalings of single coupling to get semi-couplings
	s0=GetSemiCouplingScaling(mu,_m0)
	s1=GetSemiCouplingScaling(nu,_m1)

	sMin=min(np.min(s0),np.min(s1))
	diagS0=Aux.GetCSRDiagonalMatrix(s0/sMin)
	diagS1=Aux.GetCSRDiagonalMatrix(s1/sMin)
	#diagS0=Aux.GetCSRDiagonalMatrix(s0)
	#diagS1=Aux.GetCSRDiagonalMatrix(s1)
	pi0=diagS0.dot(pi)
	pi0.sort_indices()
	pi1=pi.dot(diagS1)
	pi1.sort_indices()
	return (pi0,pi1,sMin)




def ThreshParticles(particles,minMass):
	keepPos=np.nonzero((particles[2]>minMass)+(particles[3]>minMass))
	particles=[p[keepPos] for p in particles]
	return particles
    
def GetParticles(pi,mu,nu,minMass=None):
	"""Computes list of atomic particles in WF interpolation.
	
	Arguments:
	pi: soft-marginal coupling
	mu, nu: marginals
	minMass: threshold for mass of particles.
	
	Returns:
	(particles,particlesStationary0,particlesStationary1)
	
	particles=(id0,id1,m0,m1)
	Travelling particles, id0: id of initial xpos, id1: id of final ypos, m0, m1: initial and final mass
	particlesStationary0=(id0,id1,m0,0)
	Particles that are stationary and disappear.
	id0, id1: initial and final xpos, m0: initial mass
	
	particlesStationary1=(id0,id1,0,m1)
	Analogous for second marginal.
	"""
	
	m0,m1=Common.GetMarginals(pi)
	pi0,pi1,sMin=GetSemiCouplings(pi,mu,nu,m0,m1)
	posArray=OTTools.GetFullPosArray(pi.data,pi.indices,pi.indptr)
    
	particles=(posArray[:,0],posArray[:,1],pi0.data*sMin,pi1.data*sMin)
    
	if minMass is not None:
		# threshold mass of particles if desired, to reduce computational cost of interpolation
		particles=ThreshParticles(particles,minMass)

	# determine "stationary particles" (mass that only grows or shrinks but does not move)
	m00,m01=Common.GetMarginals(pi0)
	m0Missing=mu-m00
	posArray0=np.arange(mu.shape[0])
	particlesMissing0=(posArray0,posArray0,m0Missing,np.zeros(mu.shape,dtype=np.double))
	if minMass is not None:
		particlesMissing0=ThreshParticles(particlesMissing0,minMass)
        
	m10,m11=Common.GetMarginals(pi1)
	m1Missing=nu-m11
	posArray1=np.arange(nu.shape[0])
	particlesMissing1=(posArray1,posArray1,np.zeros(nu.shape,dtype=np.double),m1Missing)
	if minMass is not None:
		particlesMissing1=ThreshParticles(particlesMissing1,minMass)
        
	return (particles,particlesMissing0,particlesMissing1)



def interpolateEuclidean(particles,xpos,ypos,t,FR_kappa):
	"""Gives Wasserstein-Fisher-Rao interpolation of a list of Dirac particles in Euclidean space.
	
	Arguments:
	particles=lists of (x0id,x1id,m0,m1), initial and final location-ids and masses
	xpos=list of positions for x0id
	ypos=list of positions for x1id
	
	Return list of (x,m):
	x: intermediate positions
	m: intermediate masses"""
	# particles is a tuple of four lists of equal length: (x0,x1,m0,m1),
	# giving initial and final location and mass of each particle
	
	x0=xpos[particles[0]]
	x1=ypos[particles[1]]
	m0=particles[2]
	m1=particles[3]
	
	if t==0:
		return (x0,m0)
	if t==1:
		return (x1,m1)
        
	# compute distances
	deltaX=np.linalg.norm(x0-x1,axis=1)

	# call distance interpolation function
	(x,m)=interpolateDistances((deltaX,m0,m1),t,FR_kappa)

	# interpolate Euclidean
	x=x0+np.einsum(x,[0],(x1-x0),[0,1],[0,1])
    
    	# return list of intermediate positions and masses at time t
	return (x,m)


def interpolateDistances(particles,t,FR_kappa):
	"""Gives Wasserstein-Fisher-Rao interpolation of a list of Dirac particles.
	
	Arguments:
	particles= lists of (dist,m0,m1)
	
	Return list of (x,m):
	x: fraction of distances travelled,
	m: intermediate masses."""

	deltaX=particles[0]	
	m0=particles[1]
	m1=particles[2]
	
	if t==0:
		return (np.full(m0.shape,0.,dtype=np.double),m0)
	if t==1:
		return (np.full(m0.shape,1.,dtype=np.double),m1)
        
	# compute aux constants
	omega=2*np.sqrt(m0*m1)*np.sin(deltaX/(2*FR_kappa)) # actually: omega / delta
	A=m0+m1-2*np.sqrt(m0*m1-omega**2/4)
	B=m0-np.sqrt(m0*m1-omega**2/4)

	# interpolate mass
	m=A*t**2-2*B*t+m0

	# interpolate location
	# initialize with start location
	x=np.full(m0.shape,0.,dtype=np.double)
	# find non-zero moment elements
	posMove=np.nonzero(omega)

	alpha=2*B[posMove]/omega[posMove]
	beta=2*(A[posMove]*t-B[posMove])/omega[posMove]
	# add movement to moving elements (multiply by normalized direction vector at the end)
	x[posMove]+=2*FR_kappa*(np.arctan(beta)+np.arctan(alpha))/deltaX[posMove]
    
    	# return list of intermediate positions and masses at time t
	return (x,m)
