import numpy as np
import scipy
import datetime

import OTTools.HierarchicalPartition as HierarchicalPartition
import Solvers.CostFunctionComputation as SolverCFC
from . import CInterface, Aux

##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
# Setup Kernel

def GetKernel_SparseCSRSanityCheck(data,indices,indptr,zeroThresh=1E-14):
	# check that each row has at least one entry
	lenList=indptr[1:]-indptr[:-1]
	if np.min(lenList)==0:
		raise ValueError("Empty row in kernel.")
	
	# check that values are not nan or inf
	countNan=np.count_nonzero(data==np.nan)
	if countNan>0:
		raise ValueError("nan in kernel.")
	
	# inf
	countInf=np.count_nonzero(data==np.inf)
	if countInf>0:
		raise ValueError("inf in kernel.")
	
	
	# check for practically-zero-rows
#	for i in range(len(indptr)-1):
#		if np.sum(data[indptr[i]:indptr[i+1]])<zeroThresh:
#			raise ValueError("Zero row in kernel.")

def GetKernel_SparseCSR(\
		partitionList,pointerListPartition,\
		method_CostFunctionProvider,\
		level, alphaList, eps,\
		cThresh=None, kThresh=None,\
		baseMeasureX=None, baseMeasureY=None,\
		verbose=False,sanityCheck=False):

	if verbose:
		time1=datetime.datetime.now()

	dim=len(alphaList)
	res=[alpha.shape[0] for alpha in alphaList]

	# hierarchical potentials
	alphaMaxList=[SolverCFC.GetSignalMinList(partitionList[i],pointerListPartition[i],alphaList[i],\
			mode=SolverCFC.MODE_MAX,lBottom=level) \
			for i in range(dim)]

	pointerAlphaList=[HierarchicalPartition.getListPointer(alphaMax) for alphaMax in alphaMaxList]
	pointerAlpha=HierarchicalPartition.getListPointer(pointerAlphaList)

	# setup cost function provider
	CFPdata=method_CostFunctionProvider(level,pointerAlpha)

	if cThresh is not None:
		slackParam=cThresh
	elif kThresh is not None:
		slackParam=-eps*np.log(kThresh)
	else:
		raise ValueError("Either cThresh oder kThresh must be provided.")
	
	# get relevant cVars
	kernelPointer,kernelSpecs=CInterface.CheckDualConstraints_Pos_generateKernel(\
			pointerListPartition, CFPdata[0][0],\
			slack=slackParam, lBottom=level)

	# cleanup of CFP
	CInterface.Delete_CostFunctionProvider(CFPdata[0][0])

	# verbose
	if verbose: print("\t\ttotal variables: ",kernelSpecs[0])

	# post processing of kernel
	CInterface.Kernel_Exponentiate(kernelPointer[0],eps)

	if verbose:
		time2=datetime.datetime.now()

	# pre-scale with base measures if desired
	if (baseMeasureX is not None) and (baseMeasureY is not None):
		baseMeasurePointer=HierarchicalPartition.getListPointer([baseMeasureX,baseMeasureY])
		CInterface.Kernel_Scale(kernelPointer[0],baseMeasurePointer)
	
	# get kernel from c interface
	kernel=CInterface.CollectPosVarList(kernelPointer[0],kernelSpecs,doDelete=True)
		
	# converting into sparse matrix formats
	kernelSparse=scipy.sparse.coo_matrix((kernel[0],(kernel[1][:,0],kernel[1][:,1])),shape=res)
	del kernel
	kernelRow=kernelSparse.tocsr()
	kernelCol=(kernelSparse.transpose()).tocsr()
	del kernelSparse

	# do sanity check if desired
	if sanityCheck:
		for kernel in [kernelRow,kernelCol]:
			GetKernel_SparseCSRSanityCheck(kernel.data,kernel.indices,kernel.indptr)
			
		
	if verbose:
		time3=datetime.datetime.now()
		print("\t\tevaluation: ",time2-time1," post-processing: ", time3-time2)
		
	return [kernelRow,kernelCol]


#def GetKernelSparseScaling_CSR(varList,alphaList,eps,\
#		level,\
#		partitionList,pointerListPartition,\
#		method_CostFunctionProvider,
#		baseMeasureX=None, baseMeasureY=None\
#		):

#	dim=len(alphaList)
#	res=[alpha.shape[0] for alpha in alphaList]

#	# setup cost function provider
#	CFPdata=method_CostFunctionProvider(level,None,alphaFinest=alphaList)

#	# get relevant cVars
#	kernel=CInterface.ReEvaluateVarList(pointerListPartition, CFPdata[0][0], level, varList)

#	# cleanup of CFP
#	CInterface.Delete_CostFunctionProvider(CFPdata[0][0])

#	# post processing of kernel
#	kernel[0][:]=np.exp(-kernel[0]/eps)
#	
#	# converting into sparse matrix formats
#	kernelSparse=scipy.sparse.coo_matrix((kernel[0],(kernel[1][:,0],kernel[1][:,1])),shape=res)
#	del kernel
#	kernelRow=kernelSparse.tocsr()
#	kernelCol=(kernelSparse.transpose()).tocsr()
#	del kernelSparse
#	
#	# pre-scale with base measures if desired
#	if (baseMeasureX is not None) and (baseMeasureY is not None):
#		diag0=Aux.GetCSRDiagonalMatrix(baseMeasureX)
#		diag1=Aux.GetCSRDiagonalMatrix(baseMeasureY)
#		kernelRow=(diag0.dot(kernelRow.dot(diag1)))
#		kernelCol=(diag1.dot(kernelCol.dot(diag0)))


#	return [kernelRow,kernelCol]

def GetKernel_DenseCSR(\
		partitionList,pointerListPartition,\
		method_CostFunctionProvider,\
		level, alphaList,eps,\
		baseMeasureX=None, baseMeasureY=None,\
		verbose=False\
		):

	if verbose:
		time1=datetime.datetime.now()
	
	res=[partition.cardLayers[level] for partition in partitionList]

	# hierarchical potentials
	if alphaList is not None:
		dim=len(alphaList)
		alphaMaxList=[SolverCFC.GetSignalMinList(partitionList[i],pointerListPartition[i],alphaList[i],\
				mode=SolverCFC.MODE_MAX,lBottom=level) \
				for i in range(dim)]

		pointerAlphaList=[HierarchicalPartition.getListPointer(alphaMax) for alphaMax in alphaMaxList]
		pointerAlpha=HierarchicalPartition.getListPointer(pointerAlphaList)
	else:
		pointerAlpha=None
	
	# setup cost function provider
	CFPdata=method_CostFunctionProvider(level,pointerAlpha)

	# get relevant cVars
	kernel=CInterface.GetDenseCosts_Pos(pointerListPartition,CFPdata[0][0],level)

	# cleanup of CFP
	CInterface.Delete_CostFunctionProvider(CFPdata[0][0])

	# post processing of kernel
	kernel[0][:]=np.exp(-kernel[0]/eps)
	
	if verbose:	
		time2=datetime.datetime.now()

	
	# converting into sparse matrix formats
	kernelSparse=scipy.sparse.coo_matrix((kernel[0],(kernel[1][:,0],kernel[1][:,1])),shape=res)
	del kernel
	kernelRow=kernelSparse.tocsr()
	kernelCol=(kernelSparse.transpose()).tocsr()
	del kernelSparse

	# pre-scale with base measures if desired
	if (baseMeasureX is not None) and (baseMeasureY is not None):
		diag0=Aux.GetCSRDiagonalMatrix(baseMeasureX)
		diag1=Aux.GetCSRDiagonalMatrix(baseMeasureY)
		kernelRow=(diag0.dot(kernelRow.dot(diag1)))
		kernelCol=(diag1.dot(kernelCol.dot(diag0)))

	if verbose:
		time3=datetime.datetime.now()
		print("\t\tevaluation: ",time2-time1," post-processing: ", time3-time2)

	return [kernelRow,kernelCol]


def GetKernel_DenseExplicitArray(c,\
		alphaList,eps,\
		baseMeasureX=None, baseMeasureY=None,\
		truncationThresh=None\
		):
	"""Compute dense kernel matrix for an explicitly given cost matrix c."""

	if alphaList is not None:
		res=c.shape
		cEff=c\
				-np.einsum(alphaList[0],[0],np.ones((res[1],),dtype=np.double),[1],[0,1])\
				-np.einsum(np.ones((res[0],),dtype=np.double),[0],alphaList[1],[1],[0,1])
	else:
		cEff=c
	

	# build dense kernel
	kernel=np.exp(-cEff/eps)
	
	# pre-scale with base measures if desired
	if (baseMeasureX is not None) and (baseMeasureY is not None):
		kernel=np.einsum(kernel,[0,1],baseMeasureX,[0],baseMeasureY,[1],[0,1])

	if truncationThresh is not None:
		pos=np.nonzero(kernel<truncationThresh)
		kernel[pos]=0

	# get transposed kernel
	kernelT=kernel.transpose().copy()


	return [kernel,kernelT]



def GetKernel_DenseArray(\
		partitionList,pointerListPartition,\
		method_CostFunctionProvider,\
		level, alphaList,eps,\
		baseMeasureX=None, baseMeasureY=None,\
		truncationThresh=None,verbose=False\
		):

	if verbose:
		time1=datetime.datetime.now()

	res=[partition.cardLayers[level] for partition in partitionList]
		
	# setup cost function provider
	CFPdata=method_CostFunctionProvider(level,pointerAlpha=None,alphaFinest=alphaList)

	# get relevant cVars
	kernelPos=CInterface.GetDenseCosts_Pos(pointerListPartition,CFPdata[0][0],level)

	# cleanup of CFP
	CInterface.Delete_CostFunctionProvider(CFPdata[0][0])

	# post processing of kernel
	kernelPos[0][:]=np.exp(-kernelPos[0]/eps)

	if verbose:	
		time2=datetime.datetime.now()

	kernel=kernelPos[0].reshape(res)
	if (baseMeasureX is not None) and (baseMeasureY is not None):
		kernel=np.einsum(kernel,[0,1],baseMeasureX,[0],baseMeasureY,[1],[0,1])
	
	if truncationThresh is not None:
		truncPos=np.nonzero(kernel<truncationThresh)
		kernel[truncPos]=0.

	kernelT=kernel.transpose().copy()

	if verbose:
		time3=datetime.datetime.now()
		print("\t\tevaluation: ",time2-time1," post-processing: ", time3-time2)

	return [kernel,kernelT]



##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
# Refinement

def RefinePotentials(partitionList, pointerListPartition, level, alphaList):
	alphaFineList=[SolverCFC.RefineSignal(partition,pointer,alpha,level) for partition,pointer,alpha\
			in zip(partitionList,pointerListPartition,alphaList)]

	return alphaFineList


def RefineKernel_CSR(partitionList, pointerListPartition, method_CostFunctionProvider,\
		level, varList, alphaFineList,\
		eps,\
		baseMeasureX=None, baseMeasureY=None,\
		verbose=False):
	
	res=[alpha.shape[0] for alpha in alphaFineList]
	CFPdata=method_CostFunctionProvider(level+1,None,alphaFinest=alphaFineList)

	# verbose
	if verbose:
		t1=datetime.datetime.now()
		

	# get relevant cVars
	kernelPointer,kernelSpecs=CInterface.RefineVarList_CSR_Pos_generateKernel(pointerListPartition, CFPdata[0][0], level, varList)

	# cleanup of CFP
	CInterface.Delete_CostFunctionProvider(CFPdata[0][0])

	# verbose
	if verbose: print("\t\ttotal variables: ",kernelSpecs[0])

	# post processing of kernel
	CInterface.Kernel_Exponentiate(kernelPointer[0],eps)

	if verbose:
		time2=datetime.datetime.now()

	# pre-scale with base measures if desired
	if (baseMeasureX is not None) and (baseMeasureY is not None):
		baseMeasurePointer=HierarchicalPartition.getListPointer([baseMeasureX,baseMeasureY])
		CInterface.Kernel_Scale(kernelPointer[0],baseMeasurePointer)
	
	# get kernel from c interface
	kernel=CInterface.CollectPosVarList(kernelPointer[0],kernelSpecs,doDelete=True)
		
	# converting into sparse matrix formats
	kernelSparse=scipy.sparse.coo_matrix((kernel[0],(kernel[1][:,0],kernel[1][:,1])),shape=res)
	del kernel
	kernelRow=kernelSparse.tocsr()
	kernelCol=(kernelSparse.transpose()).tocsr()
	del kernelSparse

	return [kernelRow,kernelCol]



##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
# PosPointer Datastructure


def GetKernel_SparsePosPointer(\
		partitionList,pointerListPartition,\
		method_CostFunctionProvider,\
		level, alphaList,eps,\
		cThresh=None, kThresh=None,\
		baseMeasure=None,\
		verbose=False):

	dim=len(alphaList)
	res=[alpha.shape[0] for alpha in alphaList]

	# hierarchical potentials
	alphaMaxList=[SolverCFC.GetSignalMinList(partitionList[i],pointerListPartition[i],alphaList[i],\
			mode=SolverCFC.MODE_MAX,lBottom=level) \
			for i in range(dim)]

	pointerAlphaList=[HierarchicalPartition.getListPointer(alphaMax) for alphaMax in alphaMaxList]
	pointerAlpha=HierarchicalPartition.getListPointer(pointerAlphaList)

	# setup cost function provider
	CFPdata=method_CostFunctionProvider(level,pointerAlpha)

	if cThresh is not None:
		slackParam=cThresh
	elif kThresh is not None:
		slackParam=-eps*np.log(kThresh)
	else:
		raise ValueError("Either cThresh oder kThresh must be provided.")

	# get relevant cVars
	kernelPointer,kernelSpecs=CInterface.CheckDualConstraints_Pos_generateKernel(\
			pointerListPartition, CFPdata[0][0],\
			slack=slackParam, lBottom=level)

	# cleanup of CFP
	CInterface.Delete_CostFunctionProvider(CFPdata[0][0])

	# verbose
	if verbose: print("\t\ttotal variables: ",kernelSpecs[0])

	# post processing of kernel
	CInterface.Kernel_Exponentiate(kernelPointer[0],eps)


	# pre-scale with base measures if desired
	if (baseMeasure is not None):
		baseMeasurePointer=HierarchicalPartition.getListPointer(baseMeasure)
		CInterface.Kernel_Scale(kernelPointer[0],baseMeasurePointer)

	return kernelPointer[0]


def RefineKernel_PosPointer(partitionList, pointerListPartition, method_CostFunctionProvider,\
		level, kernel, alphaFineList,\
		eps,\
		baseMeasure=None,\
		verbose=False):

	res=[alpha.shape[0] for alpha in alphaFineList]
	CFPdata=method_CostFunctionProvider(level+1,None,alphaFinest=alphaFineList)

	# verbose
	if verbose:
		t1=datetime.datetime.now()

	kernelPointer,kernelSpecs=CInterface.RefineVarList_Pos_Pos_generateKernel(\
			pointerListPartition, CFPdata[0][0], level, kernel)

	# cleanup of CFP
	CInterface.Delete_CostFunctionProvider(CFPdata[0][0])


	# verbose
	if verbose: print("\t\ttotal variables: ",kernelSpecs[0])

	# post processing of kernel
	CInterface.Kernel_Exponentiate(kernelPointer[0],eps)


	# pre-scale with base measures if desired
	if (baseMeasure is not None):
		baseMeasurePointer=HierarchicalPartition.getListPointer(baseMeasure)
		CInterface.Kernel_Scale(kernelPointer[0],baseMeasurePointer)

	return kernelPointer[0]

##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
# Iterations

def ErrorMarginL1(muList,marginals):
	result=0
	for ma,mu in zip(marginals,muList):
		result+=np.sum(np.abs(ma-mu))
	return result

def ErrorMarginLInf(muList,marginals):
	result=0
	for ma,mu in zip(marginals,muList):
		result=np.max([result,np.max(np.abs(ma-mu))])
	return result


def IterateToPrecision(kernel, alphaList, scalingList,pointerListScaling,muList,pointerListMu,eps,
		method_iterate, method_error,
		maxError,nInnerIterations=100,maxOuterIterations=10,verbose=False,scalingBound=None,scalingLowerBound=None,\
		residualList=None
		):
	"""Iterate until error small enough, max iterations are reached or scalings grow beyond a given bound.
	Return values: 0 -> converged, -1 -> max iterations exceeded, -2 -> scalings overflow."""


	error=np.inf
	nOuterIterations=0
	scalingOverflow=False
	
	result={"iterations":0}

	while (error>maxError) and (nOuterIterations<maxOuterIterations) and (not scalingOverflow):
	
		## run some iterations
		if verbose: t1=datetime.datetime.now()
		method_iterate(kernel, alphaList, scalingList, muList, pointerListScaling, pointerListMu, eps, nInnerIterations)
		if verbose: t2=datetime.datetime.now()
		
		## check error measure
		error=method_error(kernel, alphaList, scalingList, muList, pointerListScaling, pointerListMu, eps)
		
		## check if scalings exceed bound
		if scalingBound is not None:
			if residualList is None:
				scalingMax=max([np.max(scaling) for scaling in scalingList])
			else:
				scalingMax=max([np.max(scaling/residual) for scaling,residual in zip(scalingList,residualList)])
			if (scalingBound<scalingMax):
				scalingOverflow=True
		
		
		nOuterIterations+=1
		result["iterations"]+=nInnerIterations
		if verbose: print("\t\t",nOuterIterations," : ",error," (",(t2-t1).total_seconds(),")",flush=True)
	
	## determining proper return code
	if error<=maxError:
		if scalingLowerBound is not None:
			## do some additional checking of scalings
			if residualList is None:
				scalingMax=max([np.max(scaling) for scaling in scalingList])
			else:
				scalingMax=max([np.max(scaling/residual) for scaling,residual in zip(scalingList,residualList)])
			if scalingMax<scalingLowerBound:
				## everything went well
				result["status"]=0
			else:
				## converged but scalings are potentially large
				result["status"]=-3
		else:
			result["status"]=0
	
	elif scalingOverflow:
		## scalings too large. re-determine kernel
		result["status"]=-2
	else:
		## did not converge within specified number of max iterations
		result["status"]=-1
	
	return result

# wrapper method

def Method_IterateToPrecision(status,data,**kwargs):
	result=IterateToPrecision(\
			data["kernel"],data["alpha"],data["scaling"],data["pointerListScaling"],data["mu"],data["pointerListMu"],\
			data["eps"],\
			**kwargs)
	status["iteration_stop"]=result["status"]
	status["iterations"]+=result["iterations"]
	return (status,data)


# Absorption

def Method_AbsorbScalings(alphaList,scalingList,eps,residualScaling=None,minAlpha=None, verbose=False):

	if verbose:
		maxScaling=0
		for sl in scalingList:
			maxScaling=max(maxScaling,np.max(np.log10(sl)))
		print("\t\tmax absorbed scaling log: ",maxScaling," (eps scaled: ", eps*maxScaling," )",flush=True)

	# absorb relative to base scaling
	for i in range(len(alphaList)):
		if np.min(scalingList[i])==0:
			# if some scalings are zero, intervene
			
			if minAlpha is None:
				# throw error if no minAlpha given
				raise ValueError("Zero scaling, but no minAlpha given, at i="+str(i))
				
			nonZeroPos=np.nonzero(scalingList[i])
			zeroPos=np.nonzero(scalingList[i]==0)
			alphaList[i][zeroPos]=minAlpha[i]
			
			if residualScaling is not None:
				alphaList[i][nonZeroPos]=alphaList[i][nonZeroPos]+\
						eps*np.log(scalingList[i][nonZeroPos]/residualScaling[i][nonZeroPos])
			
				scalingList[i][nonZeroPos]=residualScaling[i][nonZeroPos]
			else:
				alphaList[i][nonZeroPos]=alphaList[i][nonZeroPos]+\
						eps*np.log(scalingList[i][nonZeroPos])
			
				scalingList[i][nonZeroPos]=1.
				
		else:
			if residualScaling is not None:
				alphaList[i][:]=alphaList[i]+eps*np.log(scalingList[i]/residualScaling[i])
				scalingList[i][:]=residualScaling[i]
			else:
				alphaList[i][:]=alphaList[i]+eps*np.log(scalingList[i])
				scalingList[i][:]=1.
								
	
	if minAlpha is not None:
		# if minAlpha constraint is imposed
		for i in range(len(alphaList)):
			alphaNew=np.maximum(minAlpha[i],alphaList[i])
			if np.sum(alphaNew>alphaList[i])>0:
				# if actually some values were truncated
				scalingList[i][:]*=np.exp((alphaList[i]-alphaNew)/eps)
				alphaList[i][:]=alphaNew


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
# Update

def Update(status,data,epsList,\
		method_getKernel, method_deleteKernel, method_absorbScaling,\
		verbose=False,absorbFinalIteration=False,maxRepeats=2):

	if status is None:
		# initial call
		if verbose: print("\tInitializing")
		status={"iterate":True, "abort" : False, "eps_n": 0, "iterations":0, "repeat":0}

		data["eps"]=epsList[0]
		status["eps"]=data["eps"]

		if data["kernel"] is None:
			# computing initial kernel
			if verbose: print("\tComputing initial kernel.")
			data["kernel"]=method_getKernel(data["kernel"],data["alpha"],data["eps"])

		return (status,data)
	
	elif status["iteration_stop"]==-3:
		# Sinkhorn converged, but scaling is not sufficiently small. do another kernel-search to be safe
		if verbose: print("\tIteration Converged. Safety absorption.")

		# Count number of safety absorptions to abort in case it does not converge.
		status["repeat"]+=1
		
		if status["repeat"]>maxRepeats:
			print("\tSafety absorptions do not seem to converge. Aborting.")
			status["iterate"]=False
			status["abort"]=True
			return (status,data)
		

		method_absorbScaling(data["alpha"],data["scaling"],data["eps"])

		# update kernel
		# delete old kernel
		method_deleteKernel(data["kernel"])
		# compute new kernel
		data["kernel"]=method_getKernel(data["kernel"],data["alpha"],data["eps"])

		return (status,data)
	
	elif status["iteration_stop"]==0:
		# Sinkhorn converged
		if status["eps_n"]<len(epsList)-1:
			# eps list not empty yet. => change epsilon
			if verbose: print("\tIteration Converged. Updating Eps.")

			# absorb old scalings
			method_absorbScaling(data["alpha"],data["scaling"],data["eps"])

			# change eps
			status["eps_n"]+=1
			data["eps"]=epsList[status["eps_n"]]
			status["eps"]=data["eps"]
			status["repeat"]=0



			# update kernel
			# delete old kernel
			method_deleteKernel(data["kernel"])
			
			# compute new kernel
			data["kernel"]=method_getKernel(data["kernel"],data["alpha"],data["eps"])
			
			return (status,data)

		else:
			# eps list empty.
			if verbose: print("\tIteration Converged. Problem solved.")

			if absorbFinalIteration:
				# absorb duals
				method_absorbScaling(data["alpha"],data["scaling"],data["eps"])

			status["iterate"]=False
			return (status,data)

	elif status["iteration_stop"]==-1:
		# Sinkhorn did not converge.
		print("\tSinkhorn didnt converge to accuarcy goal. Aborting.",flush=True)
		status["iterate"]=False
		status["abort"]=True
		return (status,data)

	elif status["iteration_stop"]==-2:
		# Sinkhorn scaling overflow. absorb duals for stabilization.
		if verbose: print("\tScaling overflow. Absorbing.")

		method_absorbScaling(data["alpha"],data["scaling"],data["eps"])

		# update kernel
		# delete old kernel
		method_deleteKernel(data["kernel"])
		# compute new kernel
		data["kernel"]=method_getKernel(data["kernel"],data["alpha"],data["eps"])

		return (status,data)

	print("\tUnknown situation occured. Aborting.")
	status["iterate"]=False
	status["abort"]=True
	return (status,data)


def SinkhornModular(method_iterate, method_update,\
		data,\
		verbose=False):

	# initialize status & data
	(status,data)=method_update(None,data)

	while status["iterate"]:
		(status,data)=method_iterate(status,data)
		(status,data)=method_update(status,data)
		if verbose: print(status,flush=True)

	return (status,data)

##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
# Complete Algorithm

#def SetupLevel(level,\
#		pre_method_getKernel, pre_method_deleteKernel, pre_method_absorbScaling, pre_method_iterate, pre_method_update,\
#		data\
#		):


#	method_getKernel = lambda kernel, scaling, alpha, eps : pre_method_getKernel(\
#		kernel, scaling, alpha, eps,\
#		level)


#	method_update = lambda status,data: pre_method_update(status,data,\
#		level,method_getKernel, pre_method_deleteKernel, pre_method_absorbScaling)

#	return ({"method_iterate": pre_method_iterate, "method_update": method_update,\
#		"data": data},\
#		{"method_getKernel" : method_getKernel,\
#				"method_deleteKernel" : pre_method_deleteKernel,\
#				"method_absorbScaling" : pre_method_absorbScaling})


def MultiscaleSolver(partitionList, pointerListPartition, muHList, epsLists,\
		get_method_getKernel, method_deleteKernel, method_absorbScaling,\
		method_iterate, get_method_update, method_refineKernel,\
		levelTop, levelBottom,\
		alphaListInit=None, kernelInit=None,\
		collectReports=False, method_getKernelVariablesCount=None,\
		refinementInterpolation=False, get_method_getKernelInterpolate=None, refinementInterpolation_steps=5,\
		verbose=False):

	level=levelTop

	
	if alphaListInit is None:
		res=np.array([partition.cardLayers[level] for partition in partitionList],dtype=np.int32)
		alphaList=[np.zeros(xres) for xres in res]
		del res
	else:
		alphaList=alphaListInit
	
	kernel=kernelInit	
	muList=[muH[level] for muH in muHList]

	if collectReports:
		statusList=[]
	
	status={"abort":False}
	totalTime=datetime.timedelta(0)
	while (level<=levelBottom) and (not status["abort"]):
	
		if verbose:
			print(level,flush=True)
		
		time1=datetime.datetime.now()

		scalingList=[np.ones_like(alpha) for alpha in alphaList]
		pointerListScaling=HierarchicalPartition.getListPointer(scalingList)
		pointerListMu=HierarchicalPartition.getListPointer(muList)

		setup={}
		setupAux={}
		setupAux["method_getKernel"]=get_method_getKernel(level, muList)
		setupAux["method_deleteKernel"]=method_deleteKernel
		setupAux["method_absorbScaling"]=method_absorbScaling

		setup["data"]={"kernel" : kernel, "alpha" : alphaList, "scaling" : scalingList, "mu" : muList,\
				"pointerListScaling" : pointerListScaling, "pointerListMu" : pointerListMu}
		setup["method_iterate"]=method_iterate
		setup["method_update"]=get_method_update(epsLists[level],\
				setupAux["method_getKernel"], setupAux["method_deleteKernel"], setupAux["method_absorbScaling"])

		if collectReports:
			variablesCountInitial=method_getKernelVariablesCount(kernel)

		(status,data)=SinkhornModular(verbose=verbose,**setup)

		if collectReports:
			variablesCountFinal=method_getKernelVariablesCount(data["kernel"])

		# refine
		if (level<levelBottom) and (not status["abort"]):
		
			# interpolation refinement
			if refinementInterpolation:

				alphaList=RefinePotentials(partitionList, pointerListPartition, level, data["alpha"])
				muList=[muH[level+1] for muH in muHList]
				
				qList=np.linspace(0,1,num=refinementInterpolation_steps,endpoint=False)
				
				for i,q in enumerate(qList):
					if verbose:
						print("\tinterpolation refinement. q={:.3f}".format(q),flush=True)
					
					scalingList=[np.ones_like(alpha) for alpha in alphaList]
					pointerListScaling=HierarchicalPartition.getListPointer(scalingList)
					pointerListMu=HierarchicalPartition.getListPointer(muList)

					setup={}
					setupAux={}
					setupAux["method_getKernel"]=get_method_getKernelInterpolate(level+1, muList, q)
					setupAux["method_deleteKernel"]=method_deleteKernel
					setupAux["method_absorbScaling"]=method_absorbScaling

					setup["data"]={"kernel" : None, "alpha" : alphaList, "scaling" : scalingList, "mu" : muList,\
							"pointerListScaling" : pointerListScaling, "pointerListMu" : pointerListMu}
					setup["method_iterate"]=method_iterate
					setup["method_update"]=get_method_update([data["eps"]],\
							setupAux["method_getKernel"], setupAux["method_deleteKernel"],\
							setupAux["method_absorbScaling"])
					
					(status,data)=SinkhornModular(verbose=True,**setup)
					
					alphaList=data["alpha"]

		
			else:
				alphaList=RefinePotentials(partitionList, pointerListPartition, level, data["alpha"])
				muList=[muH[level+1] for muH in muHList]
				kernel=method_refineKernel(level, data["kernel"], alphaList, muList,\
						data["eps"])

		time2=datetime.datetime.now()
		totalTime+=(time2-time1)
		if verbose:
			print("\tlevel time: ",time2-time1, "\ttotal time:",totalTime)
			
		if collectReports:
			status["levelTime"]=(time2-time1).total_seconds()
			status["totalTime"]=totalTime.total_seconds()
			status["vars_initial"]=variablesCountInitial
			status["vars_final"]=variablesCountFinal
			
			statusList.append(status)

		level+=1
	
	result={"data" : data, "status": status, "setup" : setup, "setupAux": setupAux}
	
	if collectReports:
		result["statusList"]=statusList
	
	return result

def GetKernelVariablesCount_CSR(kernel):
	if kernel is None:
		return 0
	else:
		return kernel[0].data.shape[0]

def GetKernelVariablesCount_PosPointer(kernel):
	if kernel is not None:
		specs=CInterface.GetSpecsPosVarList(kernel)
		return specs[0]
	else:
		return 0

def GetKernelVariablesCount_Array(kernel):
	if kernel is None:
		return 0
	else:
		shape=kernel[0].shape
		return shape[0]*shape[1]

