#ifndef ShortCutSolver_CPLEX_H_
#define ShortCutSolver_CPLEX_H_

#include<stdlib.h>
#include<PythonTypes.h>
#include"TCouplingHandler.h"
#include"TVarListHandler.h"
#include"TCPLEXNetSolver.h"
#include"Interfaces-CPLEX.h"
#include"ShortCutSolver-Tools.h"

using namespace std;

const int CH_SemiDense=0;
const int CH_Sparse=1;
const int SETUP_ERROR_CHUNKNOWN=-51;

// external functions

extern "C" {



int Setup_Solver_CPLEX(TDoubleMatrix *muX, TDoubleMatrix *muY,
		TDoubleMatrix *alpha, TDoubleMatrix *beta,
		int initializeBases,
		long int couplingHandlerAddr,
		int couplingHandlerType,
		TInteger64Matrix *Pointer);

}

#endif
