#include<stdlib.h>
#include<PythonTypes.h>
#include"TLemonSolver.h"

using namespace std;

const int Lemon_Algorithm_NS=0;
const int Lemon_Algorithm_CS=1;
const int SETUP_ERROR_LEMON_ALG_UNKNOWN=-71;

template class TLemonSolver<TCouplingHandlerDense,TAlgorithm>;
template class TLemonSolver<TCouplingHandlerSemiDense,TAlgorithm>;

template class TLemonSolver<TCouplingHandlerDense,TAlgorithmScaling>;
template class TLemonSolver<TCouplingHandlerSemiDense,TAlgorithmScaling>;


extern "C" {


int solveDense(TDoubleMatrix *c, TDoubleMatrix *muX, TDoubleMatrix *muY, TDoubleMatrix *mu,
		bool getDualVariables, TDoubleMatrix *alpha, TDoubleMatrix *beta, double cScale,
		int algorithm) {

	int xres,yres,msg;

	xres=c->dimensions[0];
	yres=c->dimensions[1];

	TCouplingHandlerDense couplingHandler(xres,yres,c->data,mu->data);

	TLemonSolverBase *solver;
	switch(algorithm) {
	case Lemon_Algorithm_NS:
		TLemonSolver<TCouplingHandlerDense,TAlgorithm> *solverNS;

		if (getDualVariables) {
			solverNS = new TLemonSolver<TCouplingHandlerDense,TAlgorithm>(&couplingHandler, muX->data, muY->data,
					alpha->data, beta->data, cScale);
		} else {
			solverNS = new TLemonSolver<TCouplingHandlerDense,TAlgorithm>(&couplingHandler, muX->data, muY->data, cScale);
		}
		solver=solverNS;
		break;
	case Lemon_Algorithm_CS:
		TLemonSolver<TCouplingHandlerDense,TAlgorithmScaling> *solverCS;

		if (getDualVariables) {
			solverCS = new TLemonSolver<TCouplingHandlerDense,TAlgorithmScaling>(&couplingHandler, muX->data, muY->data,
					alpha->data, beta->data, cScale);
		} else {
			solverCS = new TLemonSolver<TCouplingHandlerDense,TAlgorithmScaling>(&couplingHandler, muX->data, muY->data, cScale);
		}

		solver=solverCS;
		break;
	default:
		return SETUP_ERROR_LEMON_ALG_UNKNOWN;
		break;
	}

	msg=solver->setup();
	if(msg!=0) {
		return msg;
	}

	msg=solver->solve();

	delete solver;
	return msg;

}



int solveSemiDense(TDoubleMatrix *c, TDoubleMatrix *muX, TDoubleMatrix *muY,
		TDoubleMatrix *mu,
	bool getDualVariables, TDoubleMatrix *alpha, TDoubleMatrix *beta,
	TInteger32Matrix *VLXindices, TInteger32Matrix *VLXindptr,
	double cScale,
	int algorithm
	) {

	int xres,yres,msg;

	xres=c->dimensions[0];
	yres=c->dimensions[1];

	TVarListHandler *xVars;
	xVars=new TVarListHandler();
	xVars->fillFromCSRIndexList(VLXindices->data, VLXindptr->data, xres, VLXindices->dimensions[0]);

	TCouplingHandlerSemiDense couplingHandler(xres,yres,c->data,mu->data,xVars);

	TLemonSolverBase *solver;

	switch(algorithm) {
	case Lemon_Algorithm_NS:
		TLemonSolver<TCouplingHandlerSemiDense,TAlgorithm> *solverNS;

		if (getDualVariables) {
			solverNS = new TLemonSolver<TCouplingHandlerSemiDense,TAlgorithm>(&couplingHandler, muX->data, muY->data,
					alpha->data,beta->data,cScale);
		} else {
			solverNS = new TLemonSolver<TCouplingHandlerSemiDense,TAlgorithm>(&couplingHandler, muX->data, muY->data,cScale);
		}
		solver=solverNS;
		break;
	case Lemon_Algorithm_CS:
		TLemonSolver<TCouplingHandlerSemiDense,TAlgorithmScaling> *solverCS;

		if (getDualVariables) {
			solverCS = new TLemonSolver<TCouplingHandlerSemiDense,TAlgorithmScaling>(&couplingHandler, muX->data, muY->data,
					alpha->data,beta->data,cScale);
		} else {
			solverCS = new TLemonSolver<TCouplingHandlerSemiDense,TAlgorithmScaling>(&couplingHandler, muX->data, muY->data,cScale);
		}
		solver=solverCS;
		break;
	default:
		return SETUP_ERROR_LEMON_ALG_UNKNOWN;
		break;
	}


	msg=solver->setup();
	if(msg!=0) {
		return msg;
	}

	msg=solver->solve();

	// Clean up memory
	delete xVars;
	delete solver;

	return msg;
}


}
