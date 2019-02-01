#ifndef TShortCutSolver_H
#define TShortCutSolver_H

#include<cstdlib>
#include<time.h>

#include<Common/Verbose.h>
#include<Common/TVarListHandler.h>
#include<ShortCutSolver/Interfaces.h>
#include<ShortCutSolver/TShieldGenerator.h>

struct TShortCutSolverReport {
	int steps;
	double objective;
	int solved;

	// Benchmark Variables
	int t_solving;
	int t_shielding;
	int n_variables;
	int n_shielding_queries;

};


class TShortCutSolver {
public:
	static const int VCHECK_PRIMAL=0;
	static const int VCHECK_DUAL=1;
	static const int VCHECK_DUAL_PROJECT=2;

	TShortCutCouplingHandlerInterface *couplingHandler;
	TShortCutSubSolverInterfaceBase *solverInterface;
	TShieldGeneratorBase *shieldGenerator;
	int violationCheckMethod;
	TShortCutSolverReport report;

	TShortCutSolver(TShortCutCouplingHandlerInterface *_couplingHandler, TShortCutSubSolverInterfaceBase *_solverInterface,
			TShieldGeneratorBase *_shieldGenerator, int _violationCheckMethod=VCHECK_PRIMAL);
	~TShortCutSolver();

	int initialize();
	int solve();
	double getObjective();
	int step();
	int step(int steps);
	//bool dualViolationCheck(bool doProjection);
	//TVarListHandler* getSupport();
};



#endif
