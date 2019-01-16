#include"TCPLEXNetSolver.h"



TCPLEXNetSolverBase::TCPLEXNetSolverBase(double *_muX, double *_muY) {
	xres=0;
	yres=0;

	muX=_muX;
	muY=_muY;

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

	CPXenv = CPXopenCPLEX(&status);
	if (status!=0) {
		return SETUP_ERROR_ENV;
	}

	CPXnet = CPXNETcreateprob(CPXenv, &status, NULL);
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

	CPXNETaddnodes(CPXenv, CPXnet, xres, muX, NULL);
	CPXNETaddnodes(CPXenv, CPXnet, yres, muY, NULL);
	solutionStatus=false;
	return 0;
}

int TCPLEXNetSolverBase::setupArcs() {
	return MSG_NOT_IMPLEMENTED;
}

int TCPLEXNetSolverBase::setupBase(int *arcStatus, int *nodeStatus) {
	int msg;
	msg=CPXNETcopybase(CPXenv,CPXnet,arcStatus,nodeStatus);
	return msg;
}


int TCPLEXNetSolverBase::setup() {
	int msg;
	//cout << "\tsetting up env..." << endl;
	msg=setupEnvironment();
	//cout << " done" << endl;
	if(msg!=0) {
		return msg;
	}
	//cout << "\tsetting up marginals..." << endl;
	msg=setupMarginals();
	//cout << " done" << endl;
	if(msg!=0) {
		return msg;
	}
	//cout << "\tsetting up arcs..." << endl;
	msg=setupArcs();
	//cout << " done" << endl;
	if(msg!=0) {
		return msg;
	}
	return 0;
}

int TCPLEXNetSolverBase::solve() {
	return MSG_NOT_IMPLEMENTED;
}

int TCPLEXNetSolverBase::extractDuals() {
	int msg;
	double *duals;
	duals = (double*) malloc(sizeof(double) * (xres+yres));
	msg=CPXNETsolution(CPXenv, CPXnet, &msg, NULL, NULL, duals, NULL, NULL);
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
	msg=CPXNETgetbase(CPXenv,CPXnet,arcStatus,nodeStatus);
	return msg;
}


int TCPLEXNetSolverBase::deleteArcs() {
	if (!setupStatus) {
		return SETUP_ERROR_NOSETUP;
	}
	CPXNETdelarcs(CPXenv, CPXnet, 0, CPXNETgetnumarcs(CPXenv, CPXnet)-1 );

	return 0;
}

int TCPLEXNetSolverBase::cleanup() {
	if (setupStatus==false) {
		// in this case nothing has to be done and we can quit right away
		return 0;
	}

	// else, do the proper releasing
	int status;

	status=CPXNETfreeprob(CPXenv, &CPXnet);
	if(status!=0) {
		return RSLT_ERROR_FREEPROB;
	}
	CPXnet=NULL;

	status=CPXcloseCPLEX(&CPXenv);
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
