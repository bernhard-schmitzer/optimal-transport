import numpy as np
import scipy
import datetime

import OTTools.HierarchicalPartition as HierarchicalPartition
import Solvers.CostFunctionComputation as SolverCFC
from . import CInterface, Aux, Core


def GetKernel_SparseCSR(partitionList, pointerListPartition, partitionCenter, pointerPartitionCenter,\
		get_method_CostFunctionProviderPair,\
		level, alphaList, eps,\
		baseMeasureOut, baseMeasureCenter,\
		cThresh=None, kThresh=None,\
		**getKernelOptions):

	n=len(partitionList)

	kernelList=[]
	kernelListT=[]

	for i in range(n):

		alphaPair=[alphaList[i],alphaList[n+i]]
		partitionPair=[partitionList[i],partitionCenter]
		partitionPairPointer=np.array([pointerListPartition[i],pointerPartitionCenter],dtype=np.int64)

		method_CostFunctionProviderPair=get_method_CostFunctionProviderPair(i)

		result=Core.GetKernel_SparseCSR(partitionPair,partitionPairPointer,\
				method_CostFunctionProviderPair,\
				level,alphaPair,eps,\
				kThresh=kThresh, cThresh=cThresh,\
				baseMeasureX=baseMeasureOut[i], baseMeasureY=baseMeasureCenter,\
				**getKernelOptions)

		kernelList.append(result[0])
		kernelListT.append(result[1])


	return [kernelList,kernelListT]


def RefineKernel_CSR(partitionList, pointerListPartition, partitionCenter, pointerPartitionCenter,\
		get_method_CostFunctionProviderPair,\
		level, kernel, alphaList, eps, baseMeasureOut, baseMeasureCenter,\
		**refineKernelOptions):

	n=len(partitionList)

	kernelList=[]
	kernelListT=[]

	for i in range(n):

		alphaPair=[alphaList[i],alphaList[n+i]]
		partitionPair=[partitionList[i],partitionCenter]
		partitionPairPointer=np.array([pointerListPartition[i],pointerPartitionCenter],dtype=np.int64)

		method_CostFunctionProviderPair=get_method_CostFunctionProviderPair(i)

		result=Core.RefineKernel_CSR(partitionPair,partitionPairPointer,\
				method_CostFunctionProviderPair,\
				level, (kernel[0][i].indices,kernel[0][i].indptr), alphaPair,\
				eps,\
				baseMeasureX=baseMeasureOut[i], baseMeasureY=baseMeasureCenter,\
				**refineKernelOptions)

		kernelList.append(result[0])
		kernelListT.append(result[1])


	return [kernelList,kernelListT]


def GetMarginals(kernel,scalingList):
	n=len(kernel[0])
	marginalsX=[u*(k.dot(v)) for k,u,v in zip(kernel[0],scalingList[:n],scalingList[n:])]
	marginalsY=[v*(kT.dot(u)) for kT,u,v in zip(kernel[1],scalingList[:n],scalingList[n:])]
	
	return [marginalsX,marginalsY]
	
def GetKernelVariablesCount_CSR(kernel):
	if kernel is None:
		return 0
	else:
		return [Core.GetKernelVariablesCount_CSR([k]) for k in kernel[0]]


def MultiscaleSolver(partitionList, pointerListPartition,\
		partitionCenter, pointerPartitionCenter,\
		muHList, muHCenterList,\
		epsLists,\
		get_method_getKernel, method_deleteKernel, method_absorbScaling,\
		method_iterate, get_method_update, method_refineKernel,\
		levelTop, levelBottom,\
		alphaListInit=None, kernelInit=None,\
		collectReports=False, method_getKernelVariablesCount=None,\
		verbose=False):

	level=levelTop
	n=len(partitionList)


	if alphaListInit is None:
		res=np.array([partition.cardLayers[level] for partition in partitionList],dtype=np.int32)
		resCenter=partitionCenter.cardLayers[level]
		
		alphaList=[np.zeros(xres,dtype=np.double) for xres in res]\
				+[np.zeros(resCenter,dtype=np.double) for i in range(n)]

		del res
		del resCenter
		
	else:
		alphaList=alphaListInit

	kernel=kernelInit
	muList=[muH[level] for muH in muHList]
	muCenter=muHCenterList[level]

	if collectReports:
		statusList=[]

	status={"abort":False}
	totalTime=datetime.timedelta(0)
	
	while (level<=levelBottom) and (not status["abort"]):

		if verbose:
			print(level,flush=True)

		time1=datetime.datetime.now()

		scalingList=[np.ones_like(alpha) for alpha in alphaList]
		
		#pointerListScaling=HierarchicalPartition.getListPointer(scalingList)
		#pointerListMu=HierarchicalPartition.getListPointer(muList)
		pointerListScaling=None
		pointerListMu=None

		setup={}
		setupAux={}
		setupAux["method_getKernel"]=get_method_getKernel(level, muList, muCenter)
		setupAux["method_deleteKernel"]=method_deleteKernel
		setupAux["method_absorbScaling"]=method_absorbScaling

		setup["data"]={"kernel" : kernel, "alpha" : alphaList, "scaling" : scalingList, "mu" : muList,\
				"pointerListScaling" : pointerListScaling, "pointerListMu" : pointerListMu}
		setup["method_iterate"]=method_iterate
		setup["method_update"]=get_method_update(epsLists[level],\
			setupAux["method_getKernel"], setupAux["method_deleteKernel"], setupAux["method_absorbScaling"])

		if collectReports:
			variablesCountInitial=method_getKernelVariablesCount(kernel)

		(status,data)=Core.SinkhornModular(verbose=verbose,**setup)

		if collectReports:
			variablesCountFinal=method_getKernelVariablesCount(data["kernel"])

		# refine
		if (level<levelBottom) and (not status["abort"]):
			alphaList=Core.RefinePotentials(partitionList,pointerListPartition,\
					level, data["alpha"][:n])+\
					[Core.RefinePotentials([partitionCenter],[pointerPartitionCenter],\
							level, [data["alpha"][n+i]])[0] for i in range(n)]

			muList=[muH[level+1] for muH in muHList]
			muCenter=muHCenterList[level+1]

			kernel=method_refineKernel(level, data["kernel"], alphaList, muList, muCenter,\
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

