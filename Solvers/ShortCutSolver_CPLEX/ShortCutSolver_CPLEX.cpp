#include"ShortCutSolver_CPLEX.h"


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////



int Setup_Solver_CPLEX(TDoubleMatrix *muX, TDoubleMatrix *muY,
		TDoubleMatrix *alpha, TDoubleMatrix *beta,
		int initializeBases,
		long int couplingHandlerAddr,
		int couplingHandlerType,
		TInteger64Matrix *Pointer) {
	int msg;

	TCouplingHandlerExtBase *couplingHandlerInterfaceBase;
	TCPLEXNetSolverBase *solverBase;

	switch(couplingHandlerType) {
	case CH_SemiDense:
		TCouplingHandlerExt<TCouplingHandlerSemiDense> *couplingHandlerInterfaceSD;
		couplingHandlerInterfaceSD = (TCouplingHandlerExt<TCouplingHandlerSemiDense>*) couplingHandlerAddr;
		couplingHandlerInterfaceBase = couplingHandlerInterfaceSD;

		// cplex solver setup
		TCPLEXNetSolver<TCouplingHandlerSemiDense> *solverSD;
		solverSD=new TCPLEXNetSolver<TCouplingHandlerSemiDense>(couplingHandlerInterfaceSD->couplingHandler,
				muX->data,muY->data,alpha->data,beta->data);
		solverBase=solverSD;

		break;

	case CH_Sparse:
		TCouplingHandlerExt<TCouplingHandlerSparse> *couplingHandlerInterfaceSp;
		couplingHandlerInterfaceSp = (TCouplingHandlerExt<TCouplingHandlerSparse>*) couplingHandlerAddr;
		couplingHandlerInterfaceBase = couplingHandlerInterfaceSp;

		// cplex solver setup
		TCPLEXNetSolver<TCouplingHandlerSparse> *solverSp;
		solverSp=new TCPLEXNetSolver<TCouplingHandlerSparse>(couplingHandlerInterfaceSp->couplingHandler,
				muX->data,muY->data,alpha->data,beta->data);
		solverBase=solverSp;

		break;
	default:
		return SETUP_ERROR_CHUNKNOWN;
		break;
	}

	msg=solverBase->setup();
	TSolverInterfaceCPLEX *solverInterface;
	solverInterface= new TSolverInterfaceCPLEX(
			couplingHandlerInterfaceBase, solverBase,
			alpha->data, beta->data,
			initializeBases, true
			);

	if(msg!=0) {
		return msg;
	}

	// return pointers
	Pointer->data[0]=(long int) solverInterface;

	return 0;
}

