#ifndef TSinkhornSolver_H_
#define TSinkhornSolver_H_

#include<Common/PythonTypes.h>
#include<Common/ErrorCodes.h>
#include<Common/Tools.h>
#include<Common/Verbose.h>

#include<Common/THierarchicalPartition.h>
#include<Sinkhorn/TSinkhornKernel.h>

// abstract base class for sparse multi-scale sinkhorn type algorithm with eps scaling
// only implements abstract high level functions, such as:
//	iterate over scales
//	for each scale iterate over eps list
//	for each scale, eps: solve
//	during solve: check for absorption and errors

class TSinkhornSolverBase {
public:

	static constexpr int MSG_ABSORB_REITERATE=1;
	static constexpr int MSG_EXCEEDMAXITERATIONS=30101;
	static constexpr int MSG_NANSCALING=30102;	
	static constexpr int MSG_NANINERROR=30103;	
	static constexpr int MSG_ABSORB_TOOMANYABSORPTIONS=30201;

	// summarize "general" configuration parameters into one object
	struct TSinkhornSolverParameters {
		double maxError; // error accuracy to be achieved
		int maxIterations; // maximal iterations per singleSolve
		int innerIterations; // nr of iterations between absorption and error checks
		int maxAbsorptionLoops; // nr of maximal loops of absorption until algorithm must stabilize
		double absorption_scalingBound; // maximal value of a scaling variable
		double absorption_scalingLowerBound; // value above which we do a safety absorption // NOT YET IMPLEMENTED!
		double truncation_thresh; // truncation threshold in kernel sparsification
		bool refineKernel; // whether initial kernel should be generated via refinement on subsequent layers
	};

	TSinkhornSolverParameters cfg; // store all configuration parameters

	int nLayers; // total number of hierarchy layers
	int *nEpsList; // number of eps at each level
	double **epsLists; // eps lists for each level
	int layerCoarsest, layerFinest; // numbers of first and final layer to be solved

	int layer; // currently active layer
	double eps; // currently active eps
	
	bool kernelValid; // keeps track whether the currently held kernel matches the current problem setup (eps, absorption status etc.)


	TSinkhornSolverBase(
		int _nLayers,
		int *_nEpsList,
		double **_epsLists,
		int _layerCoarsest, int _layerFinest,
		TSinkhornSolverParameters _cfg
		);
		

	virtual int initialize();
	virtual int changeEps(const double newEps);
	virtual int refineDuals(__attribute__((unused)) const int newLayer) { return ERR_BASE_NOTIMPLEMENTED; };
	virtual int changeLayer(const int newLayer);
	virtual void updateParameters(TSinkhornSolverParameters newCfg);

	virtual int iterate(__attribute__((unused)) const int n) { return ERR_BASE_NOTIMPLEMENTED; };
	virtual int checkAbsorb(__attribute__((unused)) const double maxValue) { return ERR_BASE_NOTIMPLEMENTED; };
	virtual int absorb() { return ERR_BASE_NOTIMPLEMENTED; };
	virtual int getError(__attribute__((unused)) double * const result) { return ERR_BASE_NOTIMPLEMENTED; };
	virtual int generateKernel() { return ERR_BASE_NOTIMPLEMENTED; };
	
	virtual int refineKernel() { return ERR_BASE_NOTIMPLEMENTED; };
	
	
	int solveSingle();
	int solveLayer();
	int solve();
	
	
	virtual double scorePrimalUnreg() { return 0.; }; // careful: just a dummy
	
};

// class for default 2-marginal transport problem with fixed marginal constraints
class TSinkhornSolverStandard : public TSinkhornSolverBase {
public:

	// objects that are given from outside
	THierarchicalPartition *HPX, *HPY;
	double **muXH, **muYH;
	// reference measures for entropy regularization
	double **rhoXH, **rhoYH;
	THierarchicalCostFunctionProvider *costProvider;
	double **alphaH, **betaH;

	// objects the algorithm creates
	TSinkhornKernelGenerator *kernelGenerator;
	TKernelMatrix kernel, kernelT;
	int xres,yres;
	TMarginalVector u,v;

	// some info on current layer
	double *muX,*muY,*alpha,*beta;


	TSinkhornSolverStandard(
		int _nLayers,
		int *_nEpsList,
		double **_epsLists,
		int _layerCoarsest, int _layerFinest,
		TSinkhornSolverParameters _cfg,
		THierarchicalPartition *_HPX, THierarchicalPartition *_HPY,
		double **_muXH, double **_muYH,
		double **_rhoXH, double **_rhoYH,
		double **_alphaH, double **_betaH,		
		THierarchicalCostFunctionProvider *_costProvider
		);

	TSinkhornSolverStandard(
		int _nLayers,
		int *_nEpsList,
		double **_epsLists,
		int _layerCoarsest, int _layerFinest,
		TSinkhornSolverParameters _cfg,
		THierarchicalPartition *_HPX, THierarchicalPartition *_HPY,
		double **_muXH, double **_muYH,
		double **_alphaH, double **_betaH,		
		THierarchicalCostFunctionProvider *_costProvider
		) : TSinkhornSolverStandard(
				_nLayers,
				_nEpsList,
				_epsLists,
				_layerCoarsest, _layerFinest,
				_cfg,
				_HPX, _HPY,
				_muXH, _muYH,
				_muXH, _muYH,
				_alphaH, _betaH,		
				_costProvider
				) {};


	~TSinkhornSolverStandard();


	virtual int initialize();
	virtual int refineDuals(const int newLayer);
	virtual int changeLayer(const int newLayer);

	virtual int checkAbsorb(const double maxValue);
	virtual int absorb();
	virtual int generateKernel();
	virtual int refineKernel();
	

	// model specific
	virtual int iterate(const int n);
	virtual int getError(double * const result);
	virtual double scorePrimalUnreg();

	virtual double scoreTransportCost();
	
	std::vector<double> getMarginalX();
	std::vector<double> getMarginalY();
	void writeMarginalX(double *buffer);
	void writeMarginalY(double *buffer);

};


// class for 2-marginal transport problem with Kullback--Leibler soft marginal constraints
// weight of KL terms is kappa
class TSinkhornSolverKLMarginals : public TSinkhornSolverStandard {
public:

	static constexpr double DBL_INFINITY=1E100; // effective value for infinity
	double kappa;

	TSinkhornSolverKLMarginals(
		int _nLayers,
		int *_nEpsList,
		double **_epsLists,
		int _layerCoarsest, int _layerFinest,
		TSinkhornSolverParameters _cfg,
		THierarchicalPartition *_HPX, THierarchicalPartition *_HPY,
		double **_muXH, double **_muYH,
		double **_alphaH, double **_betaH,		
		THierarchicalCostFunctionProvider *_costProvider,
		double _kappa
		);

	// model specific
	virtual int iterate(const int n);
	virtual int getError(double * const result);
	virtual double scorePrimalUnreg();
	static double KL(const TMarginalVector& rho, const TMarginalVector& sigma, double sigmaThresh);
	static double KL(const TMarginalVector& rho, const TMarginalVector& sigma) {
		return KL(rho,sigma,0);
		}
	static double KLDual(const TMarginalVector& alpha, const TMarginalVector& sigma);

};

#endif
