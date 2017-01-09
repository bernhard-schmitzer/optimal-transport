#include"ShortCutSolver_Lemon.h"


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


// external functions

extern "C" {


int Setup_Solver_Lemon(TDoubleMatrix *muX, TDoubleMatrix *muY,
		TDoubleMatrix *alpha, TDoubleMatrix *beta,
		int dualOffset, double cScale, int algorithm,
		long int couplingHandlerAddr,
		int couplingHandlerType,
		TInteger64Matrix *Pointer) {
	int msg;

	TCouplingHandlerExtBase *couplingHandlerInterface;
	TLemonSolverBase *solver;

	switch(couplingHandlerType) {
	case CH_SemiDense:
		TCouplingHandlerExt<TCouplingHandlerSemiDense> *couplingHandlerInterfaceSD;
		couplingHandlerInterfaceSD = (TCouplingHandlerExt<TCouplingHandlerSemiDense>*) couplingHandlerAddr;
		couplingHandlerInterface = couplingHandlerInterfaceSD;

		// solver setup
		switch(algorithm) {
		case Lemon_Algorithm_CS:
			solver=new TLemonSolver<TCouplingHandlerSemiDense,TAlgorithmScaling>(couplingHandlerInterfaceSD->couplingHandler,
					muX->data,muY->data,alpha->data,beta->data, cScale);
			break;
		default:
			solver=new TLemonSolver<TCouplingHandlerSemiDense,TAlgorithm>(couplingHandlerInterfaceSD->couplingHandler,
					muX->data,muY->data,alpha->data,beta->data, cScale);
			break;

		}

		break;
	case CH_Sparse:
		TCouplingHandlerExt<TCouplingHandlerSparse> *couplingHandlerInterfaceSp;
		couplingHandlerInterfaceSp = (TCouplingHandlerExt<TCouplingHandlerSparse>*) couplingHandlerAddr;
		couplingHandlerInterface = couplingHandlerInterfaceSp;

		// solver setup
		switch(algorithm) {
		case Lemon_Algorithm_CS:
			solver=new TLemonSolver<TCouplingHandlerSparse,TAlgorithmScaling>(couplingHandlerInterfaceSp->couplingHandler,
					muX->data,muY->data,alpha->data,beta->data, cScale);
			break;
		default:
			solver=new TLemonSolver<TCouplingHandlerSparse,TAlgorithm>(couplingHandlerInterfaceSp->couplingHandler,
					muX->data,muY->data,alpha->data,beta->data, cScale);
			break;

		}

		break;
	default:
		return SETUP_ERROR_CHUNKNOWN;
		break;
	}

	msg=solver->setup();
	TSolverInterfaceLemon *solverInterface;
	solverInterface= new TSolverInterfaceLemon(
			couplingHandlerInterface, solver,
			alpha->data, beta->data,
			dualOffset, true
			);

	if(msg!=0) {
		return msg;
	}

	// return pointers
	Pointer->data[0]=(long int) solverInterface;

	return 0;

}



}
