#include"Interfaces-CPLEX.h"


/////////////////////////////////////////////////////////////////////////////////////////////
// TShortCutSubSolverInterfaceCPLEX Methods
/////////////////////////////////////////////////////////////////////////////////////////////


TShortCutSubSolverInterfaceCPLEX::TShortCutSubSolverInterfaceCPLEX(
		int _nLayers, double **_muXH, double **_muYH,
		TShortCutCouplingHandlerInterface *_couplingHandlerInterface,
		bool _initializeBases) {

	nLayers=_nLayers;
	muXH=_muXH;
	muYH=_muYH;
	
	couplingHandlerInterface=_couplingHandlerInterface;
	solver = NULL;

	initializeBases=_initializeBases;
	
	oldBasis=NULL;
}


TShortCutSubSolverInterfaceCPLEX::~TShortCutSubSolverInterfaceCPLEX() {
	if(solver!=NULL) {
		delete solver;
	}
	
	if(oldBasis!=NULL) {
		delete oldBasis;
	}
}


int TShortCutSubSolverInterfaceCPLEX::setLayer(
		int newLayer,
		double *_alpha, double *_beta) {

	int msg;

	if(newLayer>=nLayers) {
		return ERR_MULTISCALE_EXCEEDEDLEVELS;
	}
	
	layer=newLayer;



	if(solver!=NULL) {
		delete solver;
	}
	
	if(oldBasis!=NULL) {
		delete oldBasis;
	}
	
	alpha=_alpha;
	beta=_beta;
	
	solver=new TCPLEXNetSolver<TCouplingHandlerSparse>(couplingHandlerInterface->couplingHandler,
			muXH[layer],muYH[layer],alpha,beta);
			
	// initialize sub solver
	msg=solver->setup();
	if(msg!=0) {
		return msg;
	}
	

	return 0;
}

int TShortCutSubSolverInterfaceCPLEX::solve() {
	if(solver==NULL) { return ERR_MULTISCALE_SUBSOLVERNOTINITIALIZED; }
	int msg;
	msg=solver->solve();
	return msg;
}

int TShortCutSubSolverInterfaceCPLEX::extractBasisVarList(TVarListHandler **basis) {
	if(solver==NULL) { return ERR_MULTISCALE_SUBSOLVERNOTINITIALIZED; }

	// pointer to last varList for convenience
	TVarListHandler *xVars=couplingHandlerInterface->getXVars();
	
	// create new varListHandler to store basis
	TVarListHandler *result=new TVarListHandler();
	result->setupEmpty(solver->xres);
	// varlist signal on last varlist to extract basis info
	TVarListSignal<int> *arcStatus=new TVarListSignal<int>(xVars,0);
	// extracting
	int msg;
	msg=solver->extractBase(arcStatus->signal,NULL);
	if(msg!=0) {
		return msg;
	}
	// iterate over basis info
	TVarListHandler::TIterator it=xVars->iterationInitialize();
	while(xVars->iterate(&it)) {
		// if variable was part of basis
		if(arcStatus->signal[it.offset]==CPX_BASIC) {
			// add to basis, ignore duplicate check
			result->addToLine(it.x,it.y,false);
		}
	}
	// cleanup
	delete arcStatus;
	// write result to return variable
	*basis=result;
	
	return 0;
}

int TShortCutSubSolverInterfaceCPLEX::prepareUpdate(__attribute__((unused)) TVarListHandler *newXVars) {
	if(solver==NULL) { return ERR_MULTISCALE_SUBSOLVERNOTINITIALIZED; }

	if(initializeBases) {
		// extract basis as varList and write to oldBasis
		int msg;
		msg=extractBasisVarList(&oldBasis);
		if(msg!=0) {
			return msg;
		}
		// merge oldBasis into newXVars
		newXVars->merge(oldBasis);
	}
	return 0;
}

int TShortCutSubSolverInterfaceCPLEX::update() {
	if(solver==NULL) { return ERR_MULTISCALE_SUBSOLVERNOTINITIALIZED; }

	int msg;
	msg=solver->deleteArcs();
	if(msg!=0) {
		return msg;
	}
	msg=solver->setupArcs();
	if(msg!=0) {
		return msg;
	}

	if(initializeBases) {
		// initialize pseudo node status
		// observed from cplex: always first pseudo arc is chosen to be part of basis
		int *nodeStatus=(int*) calloc((solver->xres+solver->yres),sizeof(int));
		nodeStatus[0]=1;
		
		TVarListSignal<int> *newArcStatus=new TVarListSignal<int>(couplingHandlerInterface->getXVars(),0);
		// need this to be able to write to newArcStatus by absolute indices
		newArcStatus->computeOffsets();
		// iterate over old basis
		TVarListHandler::TIterator it=oldBasis->iterationInitialize();
		while(oldBasis->iterate(&it)) {
			// set all corresponding entries of newArcStatus to CPX_BASIC
			newArcStatus->write(it.x,it.y,CPX_BASIC);
		}
		msg=solver->setupBase(newArcStatus->signal,nodeStatus);
		free(nodeStatus);
		delete newArcStatus;
		delete oldBasis;
		oldBasis=NULL;
		if(msg!=0) {
				return msg;
		}
	}

	return 0;
}


double TShortCutSubSolverInterfaceCPLEX::getObjective() {
	if(solver==NULL) { return 0.; }
	return solver->getObjective();
}


// instantantiate solver class templates

//template class TCPLEXNetSolver<TCouplingHandlerSemiDense>;
template class TCPLEXNetSolver<TCouplingHandlerSparse>;


