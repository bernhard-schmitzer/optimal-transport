#ifndef TSinkhornKernel_H_
#define TSinkhornKernel_H_

#include<forward_list>
#include<algorithm>

#include<Common/PythonTypes.h>
#include<Common/ErrorCodes.h>
#include<Common/Tools.h>
#include<Common/Verbose.h>

#include<Common/THierarchicalPartition.h>
#include<Common/THierarchicalCostFunctionProvider.h>

#include <Eigen/SparseCore>


typedef Eigen::SparseMatrix<double,Eigen::RowMajor> TKernelMatrix;
typedef Eigen::Triplet<double> TKernelTriplet;
typedef std::vector<TKernelTriplet> TKernelTriplets;
typedef Eigen::VectorXd TMarginalVector;
typedef Eigen::VectorXi TMarginalVectorInt;
typedef Eigen::ArrayXd TMarginalArray;
typedef Eigen::ArrayXi TMarginalArrayInt;


class TSinkhornKernelGenerator {
public:

	// struct for hierarchical search for kernel entries
	struct TCandidate {
		int layer; // on which layer
		int z; // row/column index
		double v; // value
		bool operator<(const TSinkhornKernelGenerator::TCandidate& rhs) const { return v < rhs.v; };
	};

	typedef THierarchicalSearchList<TSinkhornKernelGenerator::TCandidate> TCandidateList;


	// Hierarchical Partitions
	THierarchicalPartition *HPX, *HPY;
	// marginal reference measures at finest layer
	const double * muX;
	const double* muY;
	bool useReferenceMeasures; // whether marginal reference measures have been given
	// hierarchical cost function provider	
	THierarchicalCostFunctionProvider *costProvider;
	// finest layer
	int layerBottom;
	
	// temp storage for non-zero kernel entries
	TKernelTriplets entries;
	// regularization strength
	double eps;
	// slack threshold for effective costs
	double slack;
	
	// safety mode settings
	bool useSafeMode,useFixDuals;
	
	TSinkhornKernelGenerator(
			THierarchicalPartition *_HPX, THierarchicalPartition *_HPY,
			const double * _muX, const double * _muY,
			THierarchicalCostFunctionProvider *_costProvider,
			int _layerBottom
			);
	virtual ~TSinkhornKernelGenerator() {};
	
	TKernelMatrix generate(double _eps, double _slack, bool safeMode, bool fixDuals);
	TKernelMatrix generate(double _eps, double _slack) { return generate(_eps,_slack,useSafeMode,useFixDuals); };
	void checkCell(int layer, int x, int y);
	TKernelMatrix refine(const TKernelMatrix &oldKernel, double _eps);
	
	inline void addEntry(const int x, const int y) {
		double value=costProvider->getCostEff(layerBottom, x, y);
		addEntry(x,y,value);
	}
	inline void addEntry(const int x, const int y, const double value) {
		double _value=exp(-value/eps);
		if(useReferenceMeasures) {
			_value*=muX[x]*muY[y];
		}
//		// DEBUG
//		if(std::isnan(_value)) {
//			eprintf("NAN value encountered at %d,%d\t%f\n",x,y,value);
//		}
//		// END DEBUG
		entries.push_back(TKernelTriplet(x,y,_value));							
	}
	
	std::vector<TSinkhornKernelGenerator::TCandidate> findKernelLine(int a, int mode);
	
};

//class THierarchicalKernelRow {
//public:
//	int size;
//	std::forward_list<THierarchicalKernelRowEntry> data;
//	
//	THierarchicalKernelRow(THierarchicalKernelRowEntry *newData, int n);
//	~THierarchicalKernelRow();
//	void merge(THierarchicalKernelRowEntry *newData, int n);
//};

void SinkhornAbsorbScaling(
		THierarchicalPartition *HP,
		double **alpha,
		TMarginalVector &u,
		int layer,
		double eps
		);


///////////////////////////////////////////////////////////////////////////////////

TSparseCSRContainer SinkhornKernelGetCSRData(const TKernelMatrix& kernel);
TSparsePosContainer SinkhornKernelGetPosData(const TKernelMatrix& kernel);

#endif
