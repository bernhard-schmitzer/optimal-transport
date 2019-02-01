#include"TCPLEXNetSolver.h"



TCPLEXNetSolverBase::TCPLEXNetSolverBase(double *_muX, double *_muY) {
	xres=0;
	yres=0;

	muX=_muX;
	muY=_muY;
	muYFlipped=NULL;
	muYFlippedStatus=false;

	useDualVariables=false;
	alpha=NULL;
	beta=NULL;

	objective=0.;
	nArcs=0;

	CPXenv=NULL;
	CPXnet=NULL;
	setupStatus=false;
	solutionStatus=false;


}

TCPLEXNetSolverBase::TCPLEXNetSolverBase(double *_muX, double *_muY, double *_alpha, double *_beta) {
	xres=0;
	yres=0;

	muX=_muX;
	muY=_muY;
	muYFlipped=NULL;
	muYFlippedStatus=false;

	useDualVariables=true;
	alpha=_alpha;
	beta=_beta;

	objective=0.;
	nArcs=0;

	CPXenv=NULL;
	CPXnet=NULL;
	setupStatus=false;
	solutionStatus=false;
	

}


int TCPLEXNetSolverBase::setupEnvironment() {
	if (setupStatus) {
		return SETUP_ERROR_ALREADYSET;
	}

	int status;

	CPXenv = CPXXopenCPLEX(&status);
	if (status!=0) {
		return SETUP_ERROR_ENV;
	}
	
//	double accuracy;
//	CPXXgetdblparam(CPXenv, CPX_PARAM_NETEPOPT, &accuracy);
//	eprintf("accuracy: %e\n",accuracy);
//	CPXXsetdblparam(CPXenv, CPX_PARAM_NETEPOPT, 1E-4);


	CPXnet = CPXXNETcreateprob(CPXenv, &status, NULL);
	if (status!=0) {
		return SETUP_ERROR_NET;
	}

	setupStatus=true;
	solutionStatus=false;
	return 0;
}

int TCPLEXNetSolverBase::setupMarginals() {
	if (!setupStatus) {
		return SETUP_ERROR_NOSETUP;
	}

	CPXXNETaddnodes(CPXenv, CPXnet, xres, muX, NULL);
	CPXXNETaddnodes(CPXenv, CPXnet, yres, muYFlipped, NULL);
	solutionStatus=false;
	return 0;
}

int TCPLEXNetSolverBase::setupArcs() {
	return ERR_BASE_NOTIMPLEMENTED;
}

int TCPLEXNetSolverBase::setupBase(int *arcStatus, int *nodeStatus) {
	int msg;
	msg=CPXXNETcopybase(CPXenv,CPXnet,arcStatus,nodeStatus);
	return msg;
}


int TCPLEXNetSolverBase::setup() {
	int msg;
	msg=setupEnvironment();
	if(msg!=0) {
		return msg;
	}
	msg=setupMarginals();
	if(msg!=0) {
		return msg;
	}
	msg=setupArcs();
	if(msg!=0) {
		return msg;
	}
	return 0;
}

int TCPLEXNetSolverBase::solve() {
	return ERR_BASE_NOTIMPLEMENTED;
}

int TCPLEXNetSolverBase::extractDuals() {
	int msg;
	double *duals;
	duals = (double*) malloc(sizeof(double) * (xres+yres));
	msg=CPXXNETsolution(CPXenv, CPXnet, &msg, NULL, NULL, duals, NULL, NULL);
	if(msg!=0) {
		free(duals);
		return msg;
	}
	doubleArrayCopy(duals,alpha,xres);
	doubleArrayCopy(duals+xres,beta,yres);
	doubleArrayScale(beta,-1,yres);

	free(duals);
	return 0;
}


int TCPLEXNetSolverBase::extractBase(int *arcStatus, int *nodeStatus) {
	int msg;
	msg=CPXXNETgetbase(CPXenv,CPXnet,arcStatus,nodeStatus);
	return msg;
}


int TCPLEXNetSolverBase::deleteArcs() {
	if (!setupStatus) {
		return SETUP_ERROR_NOSETUP;
	}
	CPXXNETdelarcs(CPXenv, CPXnet, 0, CPXXNETgetnumarcs(CPXenv, CPXnet)-1 );
	nArcs=0;

	return 0;
}

int TCPLEXNetSolverBase::cleanup() {
	// release memory for flipped y marginals
	if(muYFlippedStatus) {
		if(muYFlipped!=NULL) {
			free(muYFlipped);
		}
	}

	if (setupStatus==false) {
		// in this case nothing has to be done and we can quit right away
		return 0;
	}

	// else, do the proper releasing
	int status;

	status=CPXXNETfreeprob(CPXenv, &CPXnet);
	if(status!=0) {
		return RSLT_ERROR_FREEPROB;
	}
	CPXnet=NULL;

	status=CPXXcloseCPLEX(&CPXenv);
	if(status!=0) {
		return RSLT_ERROR_FREEENV;
	}
	CPXenv=NULL;

	setupStatus=false;

	return 0;
}

TCPLEXNetSolverBase::~TCPLEXNetSolverBase() {
	// call cleanup routine, just to be safe
	cleanup();
}


double TCPLEXNetSolverBase::getObjective() {
	return objective;
}


void TCPLEXNetSolverBase::setupMuYFlipped() {
	if(muYFlippedStatus) {
		return;
	}
	
	muYFlipped=(double*) malloc(sizeof(double)*yres);
	for(int y=0;y<yres;y++) {
		muYFlipped[y]=-muY[y];
	}
	muYFlippedStatus=true;
}
