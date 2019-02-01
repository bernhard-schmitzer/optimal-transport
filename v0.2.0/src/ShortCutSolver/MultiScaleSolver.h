#ifndef MultiScaleSolver_H_
#define MultiScaleSolver_H_

#include<Common/PythonTypes.h>
#include<Common/ErrorCodes.h>
#include<Common/Tools.h>
#include<Common/Verbose.h>

#include<Common/GridTools.h>
#include<Common/THierarchyBuilder.h>
#include<Common/THierarchicalPartition.h>
#include<Common/TCostFunctionProvider.h>

#include<ShortCutSolver/Interfaces.h>
#include<ShortCutSolver/TShieldGenerator.h>
#include<ShortCutSolver/TShortCutSolver.h>



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MultiScale Solver
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class TShortCutMultiScaleSolver {
public:
	// Hierarchical Partitions
	THierarchicalPartition *HPX, *HPY;
	int VCHECK;
	
	// pointers to ShortCutSolver components
	TCostFunctionProviderBase *costFunctionProvider;
	TShortCutCouplingHandlerInterface *couplingHandlerInterface;
	TShortCutSubSolverInterfaceBase *solverInterface;
	TShieldGeneratorBase *shieldGenerator;
	TShortCutSolver *ShortCutSolver;
	// whether to delete pointers after solving
	bool autoDeletePointers;
	
	// after solving, these variables store optimal primal and dual variables
	TVarListHandler *xVarsFinal;
	double *muFinal;
	double *alpha, *beta;
	double objective; // value of optimal objective


	int layerCoarsest; // index of first layer to consider for hierarchical solver

	
	TShortCutMultiScaleSolver(
			TCostFunctionProviderBase *_costFunctionProvider,
			TShortCutCouplingHandlerInterface *_couplingHandlerInterface,
			TShortCutSubSolverInterfaceBase *_solverInterface,
			TShieldGeneratorBase *_shieldGenerator,
			// Hierarchical Partitions
			THierarchicalPartition *_HPX, THierarchicalPartition *_HPY,
			// coarsest level to start with solving
			int _layerCoarsest,
			// Violation check method
			int _VCHECK
			);
			
	~TShortCutMultiScaleSolver();
	void cleanupShortCutComponents();
	int solve();

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

constexpr double MultiScaleRefineBasis_MassTolerance=1E-12;


int MultiScaleRefineBasis(THierarchicalPartition *HPX, THierarchicalPartition *HPY,
		TVarListHandler *xVarsC,
		bool *basisC, double *piC,
		double *muXF, double *muYF, TVarListHandler *xVarsF,
		int layerC,
		bool **basisFRes, double **piFRes
		);
int MultiScaleRefineBasis_NWCinCell(
		int *xList, int *yList, double *muXF, double *muYF, double *muXFSpent, double *muYFSpent,
		int xresLoc, int yresLoc,
		int *xActive, int *yActive,
		TVarListSignal<bool> *basisFSignal, TVarListSignal<double> *piFSignal,
		double m
		);

#endif
