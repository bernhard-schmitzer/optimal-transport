#include"Interfaces-Lemon.h"

TShortCutSubSolverInterfaceLemon::TShortCutSubSolverInterfaceLemon(
		int _nLayers, double **_muXH, double **_muYH,
		TShortCutCouplingHandlerInterface *_couplingHandlerInterface,
		double _measureScale, double _cScale,
		bool _dualOffset) {


	nLayers=_nLayers;
	muXH=_muXH;
	muYH=_muYH;
	couplingHandlerInterface=_couplingHandlerInterface;
	measureScale=_measureScale;
	cScale=_cScale;
	dualOffset=_dualOffset;

	solver=NULL;
	layer=0;
	
	setupDualOffset();
}


void TShortCutSubSolverInterfaceLemon::setupDualOffset() {
	if (dualOffset) {
		alphaExt=(double*) malloc(sizeof(double)*couplingHandlerInterface->getXres());
		betaExt=(double*) malloc(sizeof(double)*couplingHandlerInterface->getYres());
		// overwrite local alpha and beta
		int i;
		for(i=0;i<couplingHandlerInterface->getXres();i++) {
			alphaExt[i]=0;
		}
		for(i=0;i<couplingHandlerInterface->getYres();i++) {
			betaExt[i]=0;
		}
	}
}


TShortCutSubSolverInterfaceLemon::~TShortCutSubSolverInterfaceLemon() {
	if(solver!=NULL) {
		delete solver;
	}
	if (dualOffset) {
		free(alphaExt);
		free(betaExt);
	}
}

int TShortCutSubSolverInterfaceLemon::setLayer(
		int newLayer,
		double *_alpha, double *_beta) {

	int msg;

	if(newLayer>=nLayers) {
		return ERR_MULTISCALE_EXCEEDEDLEVELS;
	}
	
	layer=newLayer;
	alpha=_alpha;
	beta=_beta;

	if(solver!=NULL) {
		delete solver;
	}
	if (dualOffset) {
		free(alphaExt);
		free(betaExt);
		setupDualOffset();
	}


	solver=new TLemonSolver<TCouplingHandlerSparse,TAlgorithm>(couplingHandlerInterface->couplingHandler,
			muXH[layer],muYH[layer],alpha,beta, measureScale, cScale);
			
		
	// initialize sub solver
	msg=solver->setup();
	if(msg!=0) {
		return msg;
	}
	
	return 0;

	
}


int TShortCutSubSolverInterfaceLemon::solve() {
	if(solver==NULL) { return ERR_MULTISCALE_SUBSOLVERNOTINITIALIZED; }
	if (dualOffset) {
		return solver->solve(true,alphaExt,betaExt);
	} else {
		return solver->solve();
	}
}

int TShortCutSubSolverInterfaceLemon::prepareUpdate(__attribute__((unused)) TVarListHandler *newXVars) {
//	if (dualOffset) {
//		// overwrite local alpha and beta
//		int i;
//		for(i=0;i<couplingHandlerInterface->getXres();i++) {
//			alpha[i]+=alphaExt[i];
//		}
//		for(i=0;i<couplingHandlerInterface->getYres();i++) {
//			beta[i]+=betaExt[i];
//		}
//	}
	return 0;
}

int TShortCutSubSolverInterfaceLemon::update() {
	if(solver==NULL) { return ERR_MULTISCALE_SUBSOLVERNOTINITIALIZED; }
	solver->deleteArcs();
	if (dualOffset) {
		// extract local alpha and beta
		int i;
		for(i=0;i<couplingHandlerInterface->getXres();i++) {
			alphaExt[i]=alpha[i];
		}
		for(i=0;i<couplingHandlerInterface->getYres();i++) {
			betaExt[i]=beta[i];
		}
		solver->setupArcs(true,alphaExt,betaExt);
	} else {
		solver->setupArcs();
	}
	return 0;
}

double TShortCutSubSolverInterfaceLemon::getObjective() {
	if(solver==NULL) { return 0.; }
	return solver->getObjective();
}



////////////////

// instantiate template classes
//template class TLemonSolver<TCouplingHandlerSemiDense,TAlgorithm>;
template class TLemonSolver<TCouplingHandlerSparse,TAlgorithm>;

//template class TLemonSolver<TCouplingHandlerSemiDense,TAlgorithmScaling>;
//template class TLemonSolver<TCouplingHandlerSparse,TAlgorithmScaling>;


