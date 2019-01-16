from ScriptTools import ParameterType as ptype

def getParamsDefaultTransport():
	params={}
	params["setup_posScale"]=1.
	params["setup_totalMass"]=1.
	params["setup_constOffset"]=0.
	params["setup_doAbsorption"]=1
	params["setup_type_kernel"]="csr"
	params["hierarchy_lTop"]=3
	params["eps_boxScale_power"]=2.
	params["sinkhorn_maxOuter"]=10000
	params["sinkhorn_nInner"]=10
	params["sinkhorn_maxRepeats"]=100
	params["adaption_scalingBound"]=1E2
	params["adaption_scalingLowerBound"]=1E2
	params["sparsity_kThresh"]=1E-20
	params["model_FR_cMax"]=1E10
	
	return params


def getParamsDefaultBarycenter():
	params={}
	params["setup_posScale"]=1.
	params["setup_totalMass"]=1.
	params["setup_constOffset"]=0.
	params["hierarchy_lTop"]=3
	params["eps_boxScale_power"]=2.
	params["sinkhorn_maxOuter"]=10000
	params["sinkhorn_nInner"]=10
	params["sinkhorn_maxRepeats"]=100
	params["adaption_scalingBound"]=1E2
	params["adaption_scalingLowerBound"]=1E2
	params["sparsity_kThresh"]=1E-20
	params["model_FR_cMax"]=1E10

	return params



def getParamListsTransport():
	paramsListCommandLine=[\
			["setup_tag",ptype.string]\
			]

	paramsListCFGFile={\
			"setup_f1" : ptype.string,\
			"setup_f2" : ptype.string,\
			"setup_posScale" : ptype.real,\
			"setup_totalMass" : ptype.real,\
			"setup_constOffset" : ptype.real,\
			"setup_doAbsorption" : ptype.integer,\
			"setup_type_kernel" : ptype.string,\
			"hierarchy_depth" : ptype.integer,\
			"hierarchy_lTop" : ptype.integer,\
			"eps_boxScale" : ptype.real,\
			"eps_target" : ptype.real,\
			"eps_start" : ptype.real,\
			"eps_steps" : ptype.integer,\
			"eps_boxScale_power" : ptype.real,\
			"sparsity_kThresh" : ptype.real,\
			"adaption_scalingBound" : ptype.real,\
			"adaption_scalingLowerBound" : ptype.real,\
			"sinkhorn_error" : ptype.real,\
			"sinkhorn_maxOuter" : ptype.integer,\
			"sinkhorn_nInner" : ptype.integer,\
			"sinkhorn_maxRepeats" : ptype.integer,\
			"model_transportModel" : ptype.string,\
			"model_FR_kappa" : ptype.real,\
			"model_FR_cMax" : ptype.real\
			}
	return (paramsListCommandLine,paramsListCFGFile)


def getParamListsBarycenter():
	paramsListCommandLine=[\
			["setup_tag",ptype.string]\
			]

	paramsListCFGFile={\
			"setup_fileList" : ptype.stringList,\
			"setup_weightList" : ptype.realList,\
			"setup_centerRes" : ptype.integerList,\
			"setup_posScale" : ptype.real,\
			"setup_totalMass" : ptype.real,\
			"setup_constOffset" : ptype.real,\
			"setup_doAbsorption" : ptype.integer,\
			"setup_type_kernel" : ptype.string,\
			"hierarchy_depth" : ptype.integer,\
			"hierarchy_lTop" : ptype.integer,\
			"eps_boxScale" : ptype.real,\
			"eps_target" : ptype.real,\
			"eps_start" : ptype.real,\
			"eps_steps" : ptype.integer,\
			"eps_boxScale_power" : ptype.real,\
			"sparsity_kThresh" : ptype.real,\
			"adaption_scalingBound" : ptype.real,\
			"adaption_scalingLowerBound" : ptype.real,\
			"sinkhorn_error" : ptype.real,\
			"sinkhorn_maxOuter" : ptype.integer,\
			"sinkhorn_nInner" : ptype.integer,\
			"sinkhorn_maxRepeats" : ptype.integer,\
			"model_transportModel" : ptype.string,\
			"model_FR_kappa" : ptype.real,\
			"model_FR_cMax" : ptype.real\
			}
	return (paramsListCommandLine,paramsListCFGFile)


