#include"TLemonSolver.h"

TLemonSolverBase::TLemonSolverBase(double* _muX, double* _muY) {
		muX=_muX;
		muY=_muY;
		alpha=NULL;
		beta=NULL;
		useDualVariables=false;
		objective=0;
}

TLemonSolverBase::TLemonSolverBase(double* _muX, double* _muY, double* _alpha, double* _beta) {
	muX=_muX;
	muY=_muY;
	alpha=_alpha;
	beta=_beta;
	useDualVariables=true;
	objective=0;
}



