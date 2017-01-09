#include"Interfaces-Lemon.h"

template class TLemonSolver<TCouplingHandlerSemiDense,TAlgorithm>;
template class TLemonSolver<TCouplingHandlerSparse,TAlgorithm>;

template class TLemonSolver<TCouplingHandlerSemiDense,TAlgorithmScaling>;
template class TLemonSolver<TCouplingHandlerSparse,TAlgorithmScaling>;


TSolverInterfaceLemon::TSolverInterfaceLemon(
		TCouplingHandlerExtBase* _couplingHandler, TLemonSolverBase* _solver,
		bool _dualOffset, bool _deleteSolverOnDestroy) {
	couplingHandler=_couplingHandler;
	solver = _solver;

	dualOffset=_dualOffset;
	deleteSolverOnDestroy=_deleteSolverOnDestroy;

	setupDualOffset();
}

TSolverInterfaceLemon::TSolverInterfaceLemon(
		TCouplingHandlerExtBase* _couplingHandler, TLemonSolverBase* _solver,
		double* _alpha, double* _beta, bool _dualOffset, bool _deleteSolverOnDestroy) {
	alpha=_alpha;
	beta=_beta;
	couplingHandler=_couplingHandler;
	solver = _solver;

	dualOffset=_dualOffset;
	deleteSolverOnDestroy=_deleteSolverOnDestroy;

	setupDualOffset();
}

void TSolverInterfaceLemon::setupDualOffset() {
	if (dualOffset) {
		alphaExt=(double*) malloc(sizeof(double)*couplingHandler->getXres());
		betaExt=(double*) malloc(sizeof(double)*couplingHandler->getYres());
		// overwrite local alpha and beta
		int i;
		for(i=0;i<couplingHandler->getXres();i++) {
			alphaExt[i]=0;
		}
		for(i=0;i<couplingHandler->getYres();i++) {
			betaExt[i]=0;
		}
	}
}


TSolverInterfaceLemon::~TSolverInterfaceLemon() {
	if(deleteSolverOnDestroy && (solver!=NULL)) {
		delete solver;
	}
	if (dualOffset) {
		free(alphaExt);
		free(betaExt);
	}
}

int TSolverInterfaceLemon::solve() {
	return solver->solve();
}

int TSolverInterfaceLemon::prepareUpdate() {
	if (dualOffset) {
		// overwrite local alpha and beta
		int i;
		for(i=0;i<couplingHandler->getXres();i++) {
			alpha[i]+=alphaExt[i];
		}
		for(i=0;i<couplingHandler->getYres();i++) {
			beta[i]+=betaExt[i];
		}
	}
	return 0;
}

int TSolverInterfaceLemon::update() {
	solver->deleteArcs();
	if (dualOffset) {
		// extract local alpha and beta
		int i;
		for(i=0;i<couplingHandler->getXres();i++) {
			alphaExt[i]=alpha[i];
		}
		for(i=0;i<couplingHandler->getYres();i++) {
			betaExt[i]=beta[i];
		}
		solver->setupArcs(true,alphaExt,betaExt);
	} else {
		solver->setupArcs();
	}
	return 0;
}

double TSolverInterfaceLemon::getObjective() {
	return solver->getObjective(dualOffset,alphaExt,betaExt);
}
