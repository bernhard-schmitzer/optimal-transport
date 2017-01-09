import numpy as np
import scipy
from . import Common
import OTTools
from .. import Barycenter

def Iterate(kernel, alphaList, scalingList, muList,\
		eps, nInnerIterations,\
		weights, FRweight,\
		zeroHandling=False):

	prefac=2*(FRweight**2)

	n=len(muList)

	for nIt in range(nInnerIterations):
		# compute optimal beta
		if not zeroHandling:
			# without handling zero columns
			a=[kT.dot(u) for kT,u in zip(kernel[1],scalingList[:n])]
			phi=np.zeros_like(a[0])
			for i in range(n):
				phi+=weights[i]*((a[i])**(eps/(prefac+eps)))*np.exp(-alphaList[n+i]/(prefac+eps))
			phi=phi**(prefac/eps)

			for i in range(n):
				scalingList[n+i][:]=phi/(a[i]**(prefac/(prefac+eps)))*np.exp(-alphaList[n+i]/(prefac+eps))

		else:
			# with zero column handling
			a=[kT.dot(u) for kT,u in zip(kernel[1],scalingList[:n])]
		
			# identify zero and non-zero positions
			nonZeroList=[np.nonzero(a[i])[0] for i in range(n)]
		
			phi=np.zeros_like(a[0])
			for i in range(n):
				phi[nonZeroList[i]]+=weights[i]*((a[i][nonZeroList[i]])**(eps/(prefac+eps)))\
					*np.exp(-alphaList[n+i][nonZeroList[i]]/(prefac+eps))
			phi=phi**(prefac/eps)

		
		
			for i in range(n):
				scalingList[n+i][nonZeroList[i]]=phi[nonZeroList[i]]\
						/(a[i][nonZeroList[i]]**(prefac/(prefac+eps)))\
						*np.exp(-alphaList[n+i][nonZeroList[i]]/(prefac+eps))
		
			
		# u iterations as before
		if not zeroHandling:
			for i,(k,v) in enumerate(zip(kernel[0],scalingList[n:])):
				scalingList[i][:]=(muList[i]/(k.dot(v)))**(prefac/(prefac+eps))*np.exp(-alphaList[i]/(prefac+eps))
		else:
			for i,(k,v) in enumerate(zip(kernel[0],scalingList[n:])):
				conv=k.dot(v)
				posNonZero=np.nonzero(conv)[0]
				scalingList[i][posNonZero]=(muList[i][posNonZero]/(conv[posNonZero]))**(prefac/(prefac+eps))\
						*np.exp(-alphaList[i][posNonZero]/(prefac+eps))



def ScorePrimal(kernel,alphaList,scalingList,muList,eps,weights,FRweight,KLthresh=1E-10):
	prefac=2*(FRweight**2)
	n=len(weights)

	marginals=Barycenter.GetMarginals(kernel,scalingList)
	alphaEff=[alpha.copy() for alpha in alphaList]
	for i,u in enumerate(scalingList):
		nonZeroPos=np.nonzero(u)[0]
		alphaEff[i][nonZeroPos]+=eps*np.log(u[nonZeroPos])

	# barycenter
	rho=np.zeros_like(marginals[1][0])
	for i in range(n):
		rho+=weights[i]*marginals[1][i]

	scorePrimal=0.

	# F1
	for i in range(n):
		scorePrimal+=prefac*weights[i]*Common.ScoreKL(marginals[0][i],muList[i],KLthresh)
	# F2
	for i in range(n):
		scorePrimal+=prefac*weights[i]*Common.ScoreKL(marginals[1][i],rho,KLthresh)
	# G (without total kernel mass term)
	for i in range(n):
		scorePrimal+=weights[i]*(\
				np.sum(marginals[0][i]*alphaEff[i])\
				+np.sum(marginals[1][i]*alphaEff[n+i])\
				-eps*np.sum(marginals[0][i])\
				)

	return scorePrimal

def ScoreDual(kernel,alphaList,scalingList,muList,eps,weights,FRweight,betaViolationTolerance=1E-6):
	prefac=2*(FRweight**2)
	n=len(weights)
	scoreDual=0.

	masses=[np.sum(u*(k.dot(v))) for k,u,v in zip(kernel[0],scalingList[:n],scalingList[n:])]

	alphaEff=[alpha.copy() for alpha in alphaList]
	for i,u in enumerate(scalingList):
		nonZeroPos=np.nonzero(u)[0]
		alphaEff[i][nonZeroPos]+=eps*np.log(u[nonZeroPos])

	# beta constraint test
	#phi=np.zeros_like(scalingList[n])
	#for i in range(n):
	#	phi+=weights[i]*np.exp(-alphaEff[n+i]/(prefac))
	#if np.max(phi)>1+betaViolationTolerance:
	#	print("Warning! beta-constraint in dual violated.")
		#print(scalingList)
		#print(alphaList)

	for i in range(n):
		# F1*
		scoreDual-=weights[i]*prefac*Common.ScoreKLDual(-alphaEff[i]/(prefac),muList[i])
		# G* (without total kernel mass term)
		scoreDual-=eps*weights[i]*(masses[i])
	return scoreDual

def ScorePDGap(kernel, alphaList, scalingList, muList, eps, weightList, FRweight, KLthresh=1E-10):

	scorePrimal=ScorePrimal(kernel,alphaList,scalingList,muList,eps,\
			weightList, FRweight, KLthresh)

	scoreDual=ScoreDual(kernel,alphaList,scalingList,muList,eps,\
			weightList, FRweight)

	return scorePrimal-scoreDual

