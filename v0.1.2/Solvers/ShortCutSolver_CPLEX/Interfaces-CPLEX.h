#ifndef Interfaces_CPLEX_H_
#define Interfaces_CPLEX_H_

#include"Interfaces.h"
#include"TCPLEXNetSolver.h"

using namespace std;


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Short Cut Solver Interface
////////////////////////////////////////////////////////////////////////////////////////////////////////


class TSolverInterfaceCPLEX: public TSolverInterface {
public:
	TCPLEXNetSolverBase *solver;
	TCouplingHandlerExtBase *couplingHandler;
	bool deleteSolverOnDestroy;
	// whether or not to destroy the solver on destroy
	bool initializeBases;
	// initializeBases=true: optimal bases are extracted from CPLEX NET solver and
	//    re-sent to solver (on the variables that are kept) before re-solving.
	int *nodeStatus;
	TVarListSignal<int> *arcStatus;
	// nodeStatus and arcStatus store solver base information when intializeBases=true.
	//    arcStatus needs additional info, which arcs are selected
	TSolverInterfaceCPLEX(TCouplingHandlerExtBase *_couplingHandler,
			TCPLEXNetSolverBase *_solver,
			bool _initializeBases, bool _deleteSolverOnDestroy=false);
	TSolverInterfaceCPLEX(TCouplingHandlerExtBase *_couplingHandler,
			TCPLEXNetSolverBase *_solver, double *_alpha, double *_beta,
			bool _initializeBases, bool _deleteSolverOnDestroy=false);
	~TSolverInterfaceCPLEX();

	int solve();
	int prepareUpdate();
	int update();

	double getObjective();
};


#endif
