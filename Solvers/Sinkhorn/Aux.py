import numpy as np
import scipy


def SetupEpsScaling_Geometric(epsTarget,epsStart,epsSteps,verbose=False):

	if epsSteps>0:
		epsQ=(epsTarget/epsStart)**(1/epsSteps)
		if verbose: print("eps_q: ",epsQ)
		epsList=[epsStart*(epsQ**i) for i in range(epsSteps+1)]
	else:
		epsList=[epsTarget]
		epsQ=0
	
	result={"eps_q" : epsQ, "eps_list" : epsList}
	
	return result

def SetupEpsScaling_Harmonic(eps_target,eps_start,eps_steps,qMin=None,verbose=True):
	# difference of inverses of successive eps
	eps_c=(1/eps_target-1/eps_start)/eps_steps
	if verbose:
		print("eps_c: ",eps_c)
	eps_list=[1/(1/eps_start+k*eps_c) for k in range(eps_steps+1)]
    
	if qMin is not None:
		# insert intermediate eps values to soften descent at beginning
		k=0
		while k<len(eps_list)-1:
			if eps_list[k+1]/eps_list[k]<qMin:
				eps_list.insert(k+1,eps_list[k]*qMin)
			k+=1
    
	return {"eps_list" : eps_list, "eps_c" : eps_c}

def SetupEpsScaling_Scales(epsList,epsScales,levelTop=1,nOverlap=1):
	# divide eps list over scales
	epsLists=[[] for i in range(len(epsScales))]
	currentScale=levelTop
	curEps=epsScales[currentScale]
	for eps in epsList:
		while eps<epsScales[currentScale]:
			currentScale+=1
		epsLists[currentScale].append(eps)
		
	# create overlap in epsLists
	for i in range(1,len(epsScales)):
		if len(epsLists[i-1])>0:
			nOl=min(len(epsLists[i-1]),nOverlap)
			epsLists[i]=epsLists[i-1][-nOl:]+epsLists[i]
	
	result={"eps_lists" : epsLists}
	
	if len(result["eps_lists"][levelTop])==0:
		raise ValueError("levelTop has empty epsList. Try increasing hierarchy_lTop or increasing eps_boxScale.")
	
	return result


def ExtendEpsScales(epsScales,nLayer,nProbs):
	epsScales=epsScales[:nLayer]+[epsScales[nLayer] for i in range(nProbs)]+epsScales[nLayer+1:]
	return epsScales


def GetScaledSlackList(basicList,scale,offset):
	return [aList*scale+offset for aList in basicList]

def GetCSRDiagonalMatrix(u):
    n=u.shape[0]
    return scipy.sparse.csr_matrix((u,np.arange(n,dtype=np.int32),np.arange(n+1,dtype=np.int32)))

def nancheck(x,msg):
    if np.isnan(np.sum(x)):
        raise ValueError("NaN at " + msg)


def ListNonZeroAnd(dat):
	"""Returns list of indices where all elements of dat are non-zero.

	dat: list of 1d arrays of equal size."""
	# identify non-zero positions
	# start with nonZeroPos of first list
	nonZeroPos=np.nonzero(dat[0])[0]
	# then go through rest of lists and only ever keep nonZero pos
	for i in range(1,len(dat)):
		nonZeroPosNew=np.nonzero(dat[i][nonZeroPos])[0]
		nonZeroPos=nonZeroPos[nonZeroPosNew]

	return nonZeroPos

