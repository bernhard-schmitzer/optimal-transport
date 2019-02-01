#ifndef Interfaces_H_
#define Interfaces_H_

#include<cstdlib>
#include<Common/TCouplingHandler.h>
#include<Common/TVarListHandler.h>
#include<Common/Tools.h>

////////////////////////////////////////////////////////////////////////////////////////////////////////
// TShortCutSolver interacts with TCouplingHandler and the LP subsolver via two interface classes.
// TShortCutCouplingHandlerInterface is defined here
// The prototype class TShortCutSubSolverInterface is also defined
// the latter will be adapted to backend solvers such as CPLEX in other files
////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////
// TCouplingHandler interface
////////////////////////////////////////////////////////////////////////////////////////////////////////

class TShortCutCouplingHandlerInterface {
	public:
	static constexpr double MASS_TOLERANCE=1E-13;
	static constexpr double SLACK_TOLERANCE=1E-5;
	static constexpr double SLACK_INFINITY=1E12;
	TCouplingHandlerSparse *couplingHandler;

	// basic data
	int nLayers;
	int layer;
	
	// hierarchical data
	int *xresH, *yresH;
	

	TShortCutCouplingHandlerInterface(int *_xresH, int *_yresH, int _nLayers);
	~TShortCutCouplingHandlerInterface();

	int setLayer(int newLayer, TCostFunctionProviderBase *costFunctionProvider, TVarListHandler *xVars);
	
	int getXres();
	int getYres();
	TVarListHandler* getXVars();
	TVarListHandler* getSupport();
	TVarListSignal<double>* getSupportSignal();
	TSparseCSRContainer getCouplingCSR();
	TSparsePosContainer getCouplingPos();
	void updateXVars(TVarListHandler *_newXVars, bool keepXVars);
	bool dualViolationCheck(bool doProjection, double *alpha, double *beta);
};





////////////////////////////////////////////////////////////////////////////////////////////////////////
// SubSolver Interface
////////////////////////////////////////////////////////////////////////////////////////////////////////

class TShortCutSubSolverInterfaceBase {
public:
	double *alpha,*beta;
	TShortCutSubSolverInterfaceBase() {
		alpha=NULL;
		beta=NULL;
	}
	virtual ~TShortCutSubSolverInterfaceBase() {}

	// method to change layer
	virtual int setLayer(
			__attribute__((unused)) int newLayer,
			__attribute__((unused)) double *_alpha, __attribute__((unused)) double *_beta) {
		eprintf("ERROR: called setLayer() on TShortCutSubSolverInterface\n");
		return ERR_BASE_NOTIMPLEMENTED;
	}

	// aux methods for layer change
	virtual int prepareRefinement(
			__attribute__((unused)) int layerId
			) {
		return 0;
	}
	virtual int customizeRefinement(
			__attribute__((unused)) int layerId,
			__attribute__((unused)) TVarListHandler *xVars
			) {
		return 0;
	}
	

	// methods within one hierarchy layer
	virtual int solve();
	virtual int prepareUpdate(__attribute__((unused)) TVarListHandler *newXVars);
	virtual int update();
	virtual double getObjective();
};


#endif
