import numpy as np
import scipy
from .. import Aux

def ScoreKL(mu,rho,rhoThresh=0):
	if np.min(rho)<=rhoThresh:
		nonZeroPosRho=np.nonzero(rho>rhoThresh)
		zeroPosRho=np.nonzero(rho<=rhoThresh)
		
		if rhoThresh==0:
			# check if mu has support outside of support of rho
			if np.max(mu[zeroPosRho])>0:
				return np.inf

			muRes=mu[nonZeroPosRho]
			rhoRes=rho[nonZeroPosRho]
			nonZeroPosMu=np.nonzero(muRes)
			return np.sum(muRes[nonZeroPosMu]*(np.log(muRes[nonZeroPosMu]/rhoRes[nonZeroPosMu])))-\
					np.sum(mu)+np.sum(rho)
		else:
			rhoReg=np.maximum(rho,rhoThresh)
			nonZeroPosMu=np.nonzero(mu)
			return np.sum(mu[nonZeroPosMu]*(np.log(mu[nonZeroPosMu]/rhoReg[nonZeroPosMu])))-\
					np.sum(mu)+np.sum(rhoReg)
			
	else:
		if np.min(mu)==0:
			nonZeroPosMu=np.nonzero(mu)
			return np.sum(mu[nonZeroPosMu]*(np.log(mu[nonZeroPosMu]/rho[nonZeroPosMu])))-\
					np.sum(mu)+np.sum(rho)			
		else:
			return np.sum(mu*(np.log(mu/rho)))-np.sum(mu)+np.sum(rho)


def ScoreKLDual(a,b):
	return np.sum(b*(np.exp(a)-1))


def GetCouplingCSR(kernelRows,u,v,sort_indices=True):
	"""kernelRows is a 2-marginal kernel in scipy.sparse.csr_matrix format, u and v are scalings.
		Returns coupling pi in csr_matrix format."""
	diag0=Aux.GetCSRDiagonalMatrix(u)
	diag1=Aux.GetCSRDiagonalMatrix(v)
	pi=(diag0.dot(kernelRows.dot(diag1)))
	
	if sort_indices:
		pi.sort_indices()
	
	return pi

def GetMarginals(pi):
	"""Compute marginals of 2-marginal coupling pi in scipy.sparse.csr_matrix format. Deals with
		weird \"matrix\" format of scipy and returns numpy arrays instead."""
	m0=np.array(pi.sum(axis=1)).flatten()
	m1=np.array(pi.sum(axis=0)).flatten()
	return (m0,m1)
