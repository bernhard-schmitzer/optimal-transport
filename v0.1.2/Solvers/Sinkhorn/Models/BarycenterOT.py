import numpy as np
import scipy
from . import Common
from .. import Aux, Barycenter
import OTTools

def Iterate(kernel, alphaList, scalingList, muList,\
		eps, nInnerIterations,\
		weights, zeroHandling=False, setZeroInf=False):

	n=len(muList)

	for nIt in range(nInnerIterations):

		# v iteration

		if not zeroHandling:
			a=[-eps*np.log(kT.dot(u)) for kT,u in zip(kernel[1],scalingList[:n])]

			aMean=np.zeros_like(a[0])
			for i in range(n):
				aMean+=weights[i]*a[i]

			betaMean=np.zeros_like(aMean)
			for i in range(n):
				betaMean+=weights[i]*alphaList[n+i]

			for i in range(n):
				scalingList[n+i]=np.exp((a[i]-aMean-betaMean)/eps)

		else:
			a=[kT.dot(u) for kT,u in zip(kernel[1],scalingList[:n])]

			# identify non-zero positions
			nonZeroPos=Aux.ListNonZeroAnd(a)

			for i in range(n):
				a[i][nonZeroPos]=-eps*np.log(a[i][nonZeroPos])
			aMean=np.zeros_like(a[0])
			for i in range(n):
				aMean[nonZeroPos]+=weights[i]*a[i][nonZeroPos]

			betaMean=np.zeros_like(aMean)
			for i in range(n):
				betaMean[nonZeroPos]+=weights[i]*alphaList[n+i][nonZeroPos]

			for i in range(n):
				if setZeroInf:
						scalingList[n+i][:]=0.
				scalingList[n+i][nonZeroPos]=np.exp((a[i][nonZeroPos]-aMean[nonZeroPos]-betaMean[nonZeroPos])/eps)

		# u iterations as before
		for i,(k,v) in enumerate(zip(kernel[0],scalingList[n:])):
			scalingList[i]=muList[i]/(k.dot(v))

#def IterateIncreaseZero(kernel, alphaList, scalingList, muList,\
#		eps, nInnerIterations,\
#		weights, zeroHandling=False, deltaBeta=1., deltaAlpha=1.):

#	n=len(muList)

#	# v iteration

#	if not zeroHandling:
#		a=[-eps*np.log(kT.dot(u)) for kT,u in zip(kernel[1],scalingList[:n])]

#		aMean=np.zeros_like(a[0])
#		for i in range(n):
#			aMean+=weights[i]*a[i]

#		betaMean=np.zeros_like(aMean)
#		for i in range(n):
#			betaMean+=weights[i]*alphaList[n+i]

#		for i in range(n):
#			scalingList[n+i]=np.exp((a[i]-aMean-betaMean)/eps)

#	else:
#		a=[kT.dot(u) for kT,u in zip(kernel[1],scalingList[:n])]

#		# identify zero and non-zero positions
#		nonZeroList=[np.nonzero(a[i])[0] for i in range(n)]
#		zeroList=[np.nonzero(a[i]==0)[0] for i in range(n)]
#		
#				
#		# transform nonzero a to log domain
#		for i in range(n):
#			a[i][nonZeroList[i]]=-eps*np.log(a[i][nonZeroList[i]])
#				
#		# increase zero-beta
#		for i in range(n):
#			scalingList[n+i][zeroList[i]]*=np.exp(deltaBeta/eps)
#		# compute effective weights sum
#		effWeights=np.zeros_like(a[0])
#		for i in range(n):
#			effWeights[nonZeroList[i]]+=weights[i]
#		
#		# compute effective aMean
#		aMean=np.zeros_like(a[0])
#		for i in range(n):
#			aMean[nonZeroList[i]]+=weights[i]*a[i][nonZeroList[i]]
#		# compute effective betaMean (add all beta, because also zero-columns appear in zero-stabilized formula)
#		
#		betaMean=np.zeros_like(aMean)
#		for i in range(n):
#			betaMean+=weights[i]*alphaList[n+i]
#			betaMean[zeroList[i]]+=weights[i]*eps*np.log(scalingList[n+i][zeroList[i]])

#		for i in range(n):
#			scalingList[n+i][nonZeroList[i]]=np.exp((a[i][nonZeroList[i]]\
#					-((aMean[nonZeroList[i]]+betaMean[nonZeroList[i]])/effWeights[nonZeroList[i]])\
#					)/eps)
#		
#		# handle completely zero columns
#		#for i in range(n):
#		#	scalingList[n+i][totalZeroPos]*=np.exp(deltaBeta/eps)
#			

#	# u iterations as before
#	if not zeroHandling:
#		for i,(k,v) in enumerate(zip(kernel[0],scalingList[n:])):
#			scalingList[i]=muList[i]/(k.dot(v))
#	else:
#		a=[k.dot(v) for k,v in zip(kernel[0],scalingList[n:])]
#		nonZeroList=[np.nonzero(a[i])[0] for i in range(n)]
#		zeroList=[np.nonzero(a[i]==0)[0] for i in range(n)]

#		for i in range(n):
#			# increase zero rows
#			scalingList[i][zeroList[i]]*=np.exp(deltaAlpha/eps)
#			# normal iteration on other rows
#			scalingList[i][nonZeroList[i]]=muList[i][nonZeroList[i]]/a[i][nonZeroList[i]]


def ErrorMarginLInf(kernel, scalingList, muList, zeroHandling=False):

	n=len(muList)

	# outside marginal error
	marginals=[u*k.dot(v) for k,u,v in zip(kernel[0],scalingList[:n],scalingList[n:])]
	error=0
	for i in range(n):
		error=max(error,np.max(np.abs(marginals[i]-muList[i])))


	# central marginal error
	if not zeroHandling:
		marginals=[v*kT.dot(u) for (kT,u,v) in zip(kernel[1],scalingList[:n],scalingList[n:])]
		marginalsMean=np.zeros_like(marginals[0])
		for i in range(n):
			marginalsMean+=marginals[i]
		marginalsMean=marginalsMean/n

		errorCenter=0.
		for i in range(n):
			errorCenter=max(errorCenter,np.max(np.abs(marginals[i]-marginalsMean)))
	else:
		# with zero handling
		conv=[kT.dot(u) for kT,u in zip(kernel[1],scalingList[:n])]
		nonZeroPos=Aux.ListNonZeroAnd(conv)
		marginals=[v[nonZeroPos]*conv[nonZeroPos] for (conv,v) in zip(conv,scalingList[n:])]
		marginalsMean=np.zeros_like(marginals[0])
		for i in range(n):
			marginalsMean+=marginals[i]
		marginalsMean=marginalsMean/n

		errorCenter=0.
		for i in range(n):
			errorCenter=max(errorCenter,np.max(np.abs(marginals[i]-marginalsMean)))

	return max(errorCenter,error)
