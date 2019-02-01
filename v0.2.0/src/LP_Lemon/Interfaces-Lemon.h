#ifndef Interfaces_LEMON_H_
#define Interfaces_LEMON_H_

#include<cstdlib>

#include<LP_Lemon/TLemonSolver.h>
#include<ShortCutSolver/Interfaces.h>
#include<ShortCutSolver/MultiScaleSolver.h>


////////////////////////////////////////////////////////////////////////////////////////////////////////
// // Interface class for CPLEX LP solver to ShortCut algorithm
////////////////////////////////////////////////////////////////////////////////////////////////////////


class TShortCutSubSolverInterfaceLemon: public TShortCutSubSolverInterfaceBase {
public:
	TLemonSolverBase *solver;
	TShortCutCouplingHandlerInterface *couplingHandlerInterface;
	bool dualOffset;
	// dualOffset=true: effective costfunction cEff[x,y]=c[x,y]-alpha[x]-beta[y] is used for solving
	//	this might be a little bit faster.
	double *alphaExt, *betaExt;
	// keep the absolute dual variables which will be used for offsetting

	// basic data
	int nLayers;
	int layer;
	double measureScale, cScale;
	
	// hierarchical data
	double **muXH, **muYH;


	TShortCutSubSolverInterfaceLemon(
			int _nLayers, double **_muXH, double **_muYH,
			TShortCutCouplingHandlerInterface *_couplingHandlerInterface,
			double _measureScale, double _cScale,
			bool _dualOffset);

	~TShortCutSubSolverInterfaceLemon();

	// set / change layer
	int setLayer(
			int newLayer,
			double *_alpha, double *_beta);

	// on one layer
	void setupDualOffset();

	int solve();
	int prepareUpdate(__attribute__((unused)) TVarListHandler *newXVars);
	int update();

	double getObjective();
};


#endif
