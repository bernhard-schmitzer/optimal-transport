#ifndef ShortCutSolver_Lemon_H_
#define ShortCutSolver_Lemon_H_

#include<stdlib.h>
#include<PythonTypes.h>
#include"TCouplingHandler.h"
#include"TVarListHandler.h"
#include"TCostFunctionProvider-Dynamic.h"
#include"TLemonSolver.h"
#include"Interfaces-Lemon.h"
#include"ShortCutSolver-Tools.h"

using namespace std;

const int CH_SemiDense=0;
const int CH_Sparse=1;
const int SETUP_ERROR_CHUNKNOWN=-51;
const int Lemon_Algorithm_NS=0;
const int Lemon_Algorithm_CS=1;
const int SETUP_ERROR_LEMON_ALG_UNKNOWN=-71;

// external functions

extern "C" {

int Setup_Solver_Lemon(TDoubleMatrix *muX, TDoubleMatrix *muY,
		TDoubleMatrix *alpha, TDoubleMatrix *beta,
		int dualOffset, double cScale, int algorithm,
		long int couplingHandlerAddr,
		int couplingHandlerType,
		TInteger64Matrix *Pointer);

}

#endif
