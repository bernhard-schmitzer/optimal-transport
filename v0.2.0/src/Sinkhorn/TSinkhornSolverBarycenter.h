#ifndef TSinkhornSolverBarycenter_H_
#define TSinkhornSolverBarycenter_H_

#include<iostream>
using std::cout;
using std::endl;


#include<Common/PythonTypes.h>
#include<Common/ErrorCodes.h>
#include<Common/Tools.h>
#include<Common/Verbose.h>

#include<Common/THierarchicalPartition.h>
#include<Sinkhorn/TSinkhornKernel.h>
#include<Sinkhorn/TSinkhornSolver.h>

// class for default barycenter with fixed marginal constraints
class TSinkhornSolverBarycenter : public TSinkhornSolverBase {
public:

	static constexpr double DBL_ZEROCOLTOLERANCE=1E-100;
	static constexpr double DBL_MINSCALING=1E-20;
	// objects that are given from outside
	int nMarginals; // number of marginals
	double *weights; // weights for barycenter
	THierarchicalPartition **HP; // hierarchical partitions for marginals
	THierarchicalPartition *HPZ; // hierarchical partition for center
	double ***muH; // hierarchical measures for marginals
	double **muZH; // hierarchical measure for center
	THierarchicalCostFunctionProvider **costProvider; // cost providers between marginals and center
	double ***alphaH, ***betaH; // dual potentials

	// objects the algorithm creates
	TSinkhornKernelGenerator **kernelGenerator;
	std::vector<TKernelMatrix> kernel;
	std::vector<TKernelMatrix> kernelT;
	std::vector<int> res; // resolution or marginals
	int zres; // resolution of center
	std::vector<TMarginalVector> u,v; // list of scaling factors

	// some info on current layer
	std::vector<double*> mu, alpha, beta;
	double *muZ;


	TSinkhornSolverBarycenter(
		int _nLayers,
		int *_nEpsList,
		double **_epsLists,
		int _layerCoarsest, int _layerFinest,
		TSinkhornSolverParameters _cfg,
		int _nMarginals,
		double *_weights,
		THierarchicalPartition **_HP, THierarchicalPartition *_HPZ,
		double ***_muH, double **_muZH,
		double ***_alphaH, double ***_betaH,		
		THierarchicalCostFunctionProvider **_costProvider
		);

	virtual ~TSinkhornSolverBarycenter();


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
//	virtual double scorePrimalUnreg();

//	virtual double scoreTransportCost();
	
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// class for barycenter transport problem with Kullback--Leibler soft marginal constraints
// weight of KL terms is kappa
class TSinkhornSolverBarycenterKLMarginals : public TSinkhornSolverBarycenter {
public:

	static constexpr double DBL_INFINITY=1E100; // effective value for infinity
	static constexpr double KLTHRESH=1E-10; // threshold for computation of KL divergence
		// see TSinkhornSolverKLMarginals::KL for details
	double kappa;

	TSinkhornSolverBarycenterKLMarginals(
		int _nLayers,
		int *_nEpsList,
		double **_epsLists,
		int _layerCoarsest, int _layerFinest,
		TSinkhornSolverParameters _cfg,
		int _nMarginals,
		double *_weights,
		THierarchicalPartition **_HP, THierarchicalPartition *_HPZ,
		double ***_muH, double **_muZH,
		double ***_alphaH, double ***_betaH,		
		THierarchicalCostFunctionProvider **_costProvider,
		double _kappa
		);


	// model specific
	virtual int iterate(const int n);
	virtual int getError(double * const result);
	virtual int getScorePrimal(double * const result);
	virtual int getScoreDual(double * const result);
	//virtual double scorePrimalUnreg();

};


#endif
