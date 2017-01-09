from lib.header_script import *
import Solvers.Sinkhorn as Sinkhorn
import lib.header_params_Sinkhorn

# to run this script with a prepared configuration, execute the following from the terminal in the Examples/ directory:
# <python3.4 interpreter> SparseSinkhorn-01-OptimalTransport.py <config tag>
#
# possible prepared config tags are:
#
# cfg/Sinkhorn/CompareEnhancements/1
# cfg/Sinkhorn/CompareEnhancements/2
# cfg/Sinkhorn/CompareEnhancements/3
# cfg/Sinkhorn/CompareEnhancements/4
#
# cfg/Sinkhorn/ImageBenchmark/OT_256
# cfg/Sinkhorn/ImageBenchmark/WF_256
#
# cfg/Sinkhorn/Gaussians/OT_128_000-001
# cfg/Sinkhorn/Gaussians/WF_128_000-001


# verbosity: disable iteration output, may become really slow in notebooks
paramsVerbose={
		"solve_overview":True,\
		"solve_update":True,\
		"solve_kernel":True,\
		"solve_iterate":True\
		}
		
params=lib.header_params_Sinkhorn.getParamsDefaultTransport()
paramsListCommandLine,paramsListCFGFile=lib.header_params_Sinkhorn.getParamListsTransport()


# setup parameter managment
print("setup parameters",flush=True)

params=lib.header_params_Sinkhorn.getParamsDefaultTransport()
paramsListCommandLine,paramsListCFGFile=lib.header_params_Sinkhorn.getParamListsTransport()


# read setup_tag from command line argument


ScriptTools.getCommandLineParameters(params,paramsListCommandLine)


params["setup_cfgfile"]=params["setup_tag"]+".txt"

# load parameters from config file
params.update(ScriptTools.readParameters(params["setup_cfgfile"],paramsListCFGFile))

# interpreting some parameters

# totalMass regulates, whether marginals should be normalized or not.
if params["setup_totalMass"]<0:
	params["setup_totalMass"]=None
# finest level for multi-scale algorithm
params["hierarchy_lBottom"]=params["hierarchy_depth"]+1


print("final parameter settings")
for k in sorted(params.keys()):
	print("\t",k,params[k])
	

# define problem: setup marginals
print("problem setup",flush=True)

def loadProblem(filename):
	img=sciio.loadmat(filename)["a"]
	return img

def setupDensity(img,posScale,totalMass,constOffset,keepZero):
	(mu,pos)=OTTools.processDensity_Grid(img,totalMass=totalMass,constOffset=constOffset,keepZero=keepZero)
	pos=pos/posScale
	return (mu,pos,img.shape)

problemData=[setupDensity(loadProblem(filename),posScale=params["setup_posScale"],\
		totalMass=params["setup_totalMass"],constOffset=params["setup_constOffset"],keepZero=False)\
		for filename in [params["setup_f1"],params["setup_f2"]]]


nProblems=len(problemData)

# set up eps-scaling

# geometric scaling from eps_start to eps_target in eps_steps+1 steps
params.update(
		Sinkhorn.Aux.SetupEpsScaling_Geometric(params["eps_target"],params["eps_start"],params["eps_steps"],\
		verbose=True))

# determine finest epsilon for each hierarchy level.
# on coarsest level it is given by params["eps_boxScale"]**params["eps_boxScale_power"]
# with each level, the finest scale params["eps_boxScale"] is effectively divided by 2
params["eps_scales"]=[(params["eps_boxScale"]/(2**n))**params["eps_boxScale_power"]\
		for n in range(params["hierarchy_depth"]+1)]+[0]
print("eps_scales:\t",params["eps_scales"])

# divide eps_list into eps_lists, one for each hierarchy scale, divisions determined by eps_scales.
params.update(Sinkhorn.Aux.SetupEpsScaling_Scales(params["eps_list"],params["eps_scales"],\
		levelTop=params["hierarchy_lTop"], nOverlap=1))


## setup hierarchical partitions
print("setup hierarchical partitions",flush=True)
partitionChildMode=HierarchicalPartition.THPMode_Tree

# constructing basic partitions
partitionList=[]
for i in range(nProblems):
	partition=HierarchicalPartition.GetPartition(problemData[i][1],params["hierarchy_depth"],partitionChildMode,\
			box=None, signal_pos=True, signal_radii=True,clib=SolverCFC, export=False, verbose=False,\
			finestDimsWarning=False)
	partitionList.append(partition)

# exporting partitions
pointerListPartition=np.zeros((nProblems),dtype=np.int64)
for i in range(nProblems):
	pointerListPartition[i]=SolverCFC.Export(partitionList[i])

muHList=[SolverCFC.GetSignalMass(pointer,partition,aprob[0])
		for pointer,partition,aprob in zip(pointerListPartition,partitionList,problemData)]

# pointer lists
pointerPosList=[HierarchicalPartition.getSignalPointer(partition,"pos") for partition in partitionList]
pointerRadiiList=[HierarchicalPartition.getSignalPointer(partition,"radii",lBottom=partition.nlayers-2)
		for partition in partitionList]
pointerListPos=np.array([pointerPos.ctypes.data for pointerPos in pointerPosList],dtype=np.int64)
pointerListRadii=np.array([pointerRadii.ctypes.data for pointerRadii in pointerRadiiList],dtype=np.int64)

# print a few stats on the created problem
for i,partition in enumerate(partitionList):
	print("cells in partition {:d}: ".format(i),partition.cardLayers)
	
	
print("setup solver")
# model specific stuff
import Solvers.Sinkhorn.Models.OT as ModelOT

if params["model_transportModel"]=="ot":

	method_CostFunctionProvider = lambda level, pointerAlpha, alphaFinest=None :\
			Sinkhorn.CInterface.Setup_CostFunctionProvider_SquaredEuclidean(pointerListPos,\
					partitionList[0].ndim,level,pointerListRadii,pointerAlpha,alphaFinest)

	method_iterate_iterate = lambda kernel, alphaList, scalingList, muList, pointerListScaling, pointerListMu,\
			eps, nInnerIterations: \
			ModelOT.Iterate(kernel[0],kernel[1],scalingList[0],scalingList[1],muList[0],muList[1],nInnerIterations)

	def method_iterate_error(kernel, alphaList, scalingList, muList, pointerListScaling, pointerListMu, eps):
			return ModelOT.ErrorMarginLInf(kernel[0],kernel[1],scalingList[0],scalingList[1],muList[0],muList[1])

elif params["model_transportModel"]=="wf":

	import Solvers.Sinkhorn.Models.WF as ModelWF

	method_CostFunctionProvider = lambda level, pointerAlpha, alphaFinest=None :\
			Sinkhorn.CInterface.Setup_CostFunctionProvider_SquaredEuclideanWF(pointerListPos,\
					partitionList[0].ndim,level,pointerListRadii,pointerAlpha,alphaFinest,\
					FR_kappa=params["model_FR_kappa"],FR_cMax=params["model_FR_cMax"])

	method_iterate_iterate = lambda kernel, alphaList, scalingList, muList, pointerListScaling, pointerListMu,\
			eps, nInnerIterations: \
			ModelWF.Iterate(kernel[0],kernel[1],alphaList[0],alphaList[1],scalingList[0],scalingList[1],\
					muList[0],muList[1],eps,params["model_FR_kappa"],nInnerIterations)

	method_iterate_error = lambda kernel, alphaList, scalingList, muList, pointerListScaling, pointerListMu, eps: \
			ModelWF.ScorePDGap(kernel[0],kernel[1],alphaList[0],alphaList[1],\
					scalingList[0],scalingList[1],muList[0],muList[1],\
					eps,params["model_FR_kappa"])

else:
	raise ValueError("model_transportModel not recognized: "+params["model_transportModel"])


# data structure choice for kernel
if params["setup_type_kernel"]=="csr":
	get_method_getKernel = lambda level, muList:\
			lambda kernel, alpha, eps:\
					Sinkhorn.GetKernel_SparseCSR(
							partitionList,pointerListPartition,\
							method_CostFunctionProvider,\
							level, alpha, eps,\
							kThresh=params["sparsity_kThresh"],\
							baseMeasureX=muList[0], baseMeasureY=muList[1],\
							sanityCheck=False,\
							verbose=paramsVerbose["solve_kernel"])

	method_deleteKernel = lambda kernel : None

	method_refineKernel = lambda level, kernel, alphaList, muList, eps:\
					Sinkhorn.RefineKernel_CSR(partitionList, pointerListPartition,\
							method_CostFunctionProvider,\
							level, (kernel[0].indices,kernel[0].indptr), alphaList,\
							eps,\
							baseMeasureX=muList[0], baseMeasureY=muList[1],\
							verbose=paramsVerbose["solve_kernel"])

	method_getKernelVariablesCount=Sinkhorn.GetKernelVariablesCount_CSR

elif params["setup_type_kernel"]=="dense":
	get_method_getKernel = lambda level, muList:\
			lambda kernel, alpha, eps:\
					Sinkhorn.GetKernel_DenseArray(
							partitionList,pointerListPartition,\
							method_CostFunctionProvider,\
							level, alpha, eps,\
							baseMeasureX=muList[0], baseMeasureY=muList[1],\
							truncationThresh=1E-200,verbose=paramsVerbose["solve_kernel"]\
							)


	method_deleteKernel = lambda kernel : None

	method_refineKernel = lambda level, kernel, alphaList, muList, eps:\
			None
	
	method_getKernelVariablesCount=Sinkhorn.GetKernelVariablesCount_Array

else:
	raise ValueError("setup_type_kernel not recognized: "+params["setup_type_kernel"])



if params["setup_doAbsorption"]==1:
	method_absorbScaling = lambda alphaList,scalingList,eps:\
					Sinkhorn.Method_AbsorbScalings(alphaList,scalingList,eps,\
							residualScaling=None,minAlpha=None,verbose=False)
else:
	method_absorbScaling = lambda alphaList,scalingList,eps: None

get_method_update = lambda epsList, method_getKernel, method_deleteKernel, method_absorbScaling:\
		lambda status, data:\
				Sinkhorn.Update(status,data,epsList,\
						method_getKernel, method_deleteKernel, method_absorbScaling,\
						absorbFinalIteration=True,maxRepeats=params["sinkhorn_maxRepeats"],\
						verbose=paramsVerbose["solve_update"]\
						)


method_iterate = lambda status,data : Sinkhorn.Method_IterateToPrecision(status,data,\
				method_iterate=method_iterate_iterate,method_error=method_iterate_error,\
				maxError=params["sinkhorn_error"],\
				nInnerIterations=params["sinkhorn_nInner"],maxOuterIterations=params["sinkhorn_maxOuter"],\
				scalingBound=params["adaption_scalingBound"],scalingLowerBound=params["adaption_scalingLowerBound"],\
				verbose=paramsVerbose["solve_iterate"])
				
print("solving",flush=True)


result=Sinkhorn.MultiscaleSolver(partitionList,pointerListPartition,muHList,params["eps_lists"],\
		get_method_getKernel,method_deleteKernel,method_absorbScaling,\
		method_iterate,get_method_update,method_refineKernel,\
		params["hierarchy_lTop"],params["hierarchy_lBottom"],\
		collectReports=True,method_getKernelVariablesCount=method_getKernelVariablesCount,\
		verbose=paramsVerbose["solve_overview"],\
		)
