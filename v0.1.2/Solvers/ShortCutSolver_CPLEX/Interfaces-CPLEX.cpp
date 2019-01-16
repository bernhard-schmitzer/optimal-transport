#include"Interfaces-CPLEX.h"

// instantantiate solver class templates

template class TCPLEXNetSolver<TCouplingHandlerSemiDense>;
template class TCPLEXNetSolver<TCouplingHandlerSparse>;


/////////////////////////////////////////////////////////////////////////////////////////////
// TSolverInterfaceCPLEX Methods
/////////////////////////////////////////////////////////////////////////////////////////////


TSolverInterfaceCPLEX::TSolverInterfaceCPLEX(TCouplingHandlerExtBase *_couplingHandler,
		TCPLEXNetSolverBase *_solver,
		bool _initializeBases, bool _deleteSolverOnDestroy) {
	couplingHandler=_couplingHandler;
	solver = _solver;

	initializeBases=_initializeBases;
	deleteSolverOnDestroy=_deleteSolverOnDestroy;
	nodeStatus=NULL;
	arcStatus=NULL;
}

TSolverInterfaceCPLEX::TSolverInterfaceCPLEX(TCouplingHandlerExtBase *_couplingHandler,
		TCPLEXNetSolverBase *_solver,
		 double *_alpha, double *_beta,
			bool _initializeBases, bool _deleteSolverOnDestroy) {
	alpha=_alpha;
	beta=_beta;
	couplingHandler=_couplingHandler;
	solver = _solver;

	initializeBases=_initializeBases;
	deleteSolverOnDestroy=_deleteSolverOnDestroy;
	nodeStatus=NULL;
	arcStatus=NULL;
}

TSolverInterfaceCPLEX::~TSolverInterfaceCPLEX() {
	if(deleteSolverOnDestroy && (solver!=NULL)) {
		delete solver;
	}
}

int TSolverInterfaceCPLEX::solve() {
	int msg;
	msg=solver->solve();
	return msg;
}

int TSolverInterfaceCPLEX::prepareUpdate() {
	if(initializeBases) {
		arcStatus=new TVarListSignal<int>(new TVarListHandler(couplingHandler->getXVars()),0);
		nodeStatus=(int*) malloc(sizeof(int)*(couplingHandler->getXres()+couplingHandler->getYres()));
		int msg;
		msg=solver->extractBase(arcStatus->signal,nodeStatus);
		return msg;
	}
	return 0;
}

int TSolverInterfaceCPLEX::update() {
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
		TVarListSignal<int> *newArcStatus=new TVarListSignal<int>(couplingHandler->getXVars(),0);
		newArcStatus->transcribeSorted(arcStatus,0);
		msg=solver->setupBase(newArcStatus->signal,nodeStatus);
		free(nodeStatus);
		delete arcStatus->varList;
		delete arcStatus;
		delete newArcStatus;
		if(msg!=0) {
				return msg;
		}
	}

	return 0;
}

double TSolverInterfaceCPLEX::getObjective() {
	return solver->getObjective();
}


