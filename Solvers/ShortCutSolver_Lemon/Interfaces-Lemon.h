#ifndef Interfaces_LEMON_H_
#define Interfaces_LEMON_H_

#include"Interfaces.h"
#include"TLemonSolver.h"

using namespace std;


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Short Cut Solver Interface
////////////////////////////////////////////////////////////////////////////////////////////////////////


class TSolverInterfaceLemon: public TSolverInterface {
public:
	TLemonSolverBase *solver;
	TCouplingHandlerExtBase *couplingHandler;
	bool deleteSolverOnDestroy;
	// whether or not to destroy the solver on destroy
	bool dualOffset;
	// dualOffset=true: effective costfunction cEff[x,y]=c[x,y]-alpha[x]-beta[y] is used for solving
	//	this might be a little bit faster.
	double *alphaExt, *betaExt;
	// keep the absolute dual variables which will be used for offsetting
	TSolverInterfaceLemon(TCouplingHandlerExtBase *_couplingHandler,
			TLemonSolverBase *_solver,
			bool _dualOffset, bool _deleteSolverOnDestroy=false);
	TSolverInterfaceLemon(TCouplingHandlerExtBase *_couplingHandler,
			TLemonSolverBase *_solver, double *_alpha, double *_beta,
			bool _dualOffset, bool _deleteSolverOnDestroy=false);
	~TSolverInterfaceLemon();

	void setupDualOffset();

	int solve();
	int prepareUpdate();
	int update();

	double getObjective();
};


#endif
