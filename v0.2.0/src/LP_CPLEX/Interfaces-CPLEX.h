#ifndef Interfaces_CPLEX_H_
#define Interfaces_CPLEX_H_

#include<cstdlib>

#include<LP_CPLEX/TCPLEXNetSolver.h>
#include<ShortCutSolver/Interfaces.h>
#include<ShortCutSolver/MultiScaleSolver.h>


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Interface class for CPLEX LP solver to ShortCut algorithm
////////////////////////////////////////////////////////////////////////////////////////////////////////

class TShortCutSubSolverInterfaceCPLEX: public TShortCutSubSolverInterfaceBase {
public:
	TCPLEXNetSolverBase *solver;
	TShortCutCouplingHandlerInterface *couplingHandlerInterface;
	TVarListHandler *oldBasis;

	// basic data
	bool initializeBases;
			// initializeBases=true: optimal bases are extracted from CPLEX NET solver and
			//    re-sent to solver (on the variables that are kept) before re-solving.
	int nLayers;
	int layer;
	
	// hierarchical data
	double **muXH, **muYH;


	
	
	TShortCutSubSolverInterfaceCPLEX(
			int _nLayers, double **_muXH, double **_muYH,
			TShortCutCouplingHandlerInterface *_couplingHandlerInterface,
			bool _initializeBases);
	~TShortCutSubSolverInterfaceCPLEX();

	int setLayer(
			int newLayer,
			double *_alpha, double *_beta);


	int solve();
	int prepareUpdate(__attribute__((unused)) TVarListHandler *newXVars);
	int update();

	double getObjective();

	int extractBasisVarList(TVarListHandler **basis);

};


#endif
