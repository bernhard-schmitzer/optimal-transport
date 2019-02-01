#ifndef TShieldGenerator_H
#define TShieldGenerator_H

#include<cstdlib>
#include<cmath>

#include<Common/GridTools.h>
#include<Common/TVarListHandler.h>
#include<Common/THierarchicalPartition.h>


class TShieldGeneratorBase {
public:
	static constexpr double shieldingTolerance=1E-5;
	// benchmark variable to be used by some shielding generators
	int n_shielding_queries;
	TShieldGeneratorBase() { n_shielding_queries=0; };
	virtual ~TShieldGeneratorBase() {};
	virtual void generateShield(__attribute__((unused)) TVarListHandler *xVars,
			__attribute__((unused)) TVarListHandler *xSupport) {
		eprintf("ERROR: called generateShield() on TShieldGeneratorBase\n");
	}
	static void getXMap(int *xMap, TVarListHandler *xSupport);
	virtual int setLayer(__attribute__((unused)) int newLayer) {
		eprintf("ERROR: called setLayer() on TShieldGeneratorBase\n");
		return ERR_BASE_NOTIMPLEMENTED;
	}
};


class TShieldGeneratorGrid_SqrEuclidean : public TShieldGeneratorBase {
public:
	// basic info
	int dim; // number of spatial dimensions
	int nLayers; // number of layers
	int layer; // number of currently active layer
	
	// hierarchical data
	int *xDimH,*yDimH; // number of points along each axis
	
	// data at current layer
	int *xDims,*yDims; // xDims,yDims: number of points along each axis
	int *xStrides, *yStrides; // strides: difference in index, for moving on grid point in a given dimension



	TShieldGeneratorGrid_SqrEuclidean(int _dim, int *_xDimH, int *_yDimH, int _nLayers);
	~TShieldGeneratorGrid_SqrEuclidean();

	int setLayer(int newLayer);

	void generateShield(TVarListHandler *xVars, TVarListHandler *xSupport);
	void iterateXVariables(TVarListHandler *xVars, int *xMap, int *xPos, int d);

	void addVariables_Shields(TVarListHandler *xVars, int *xMap, int *xPos);
	void addVariables_Rectangles(TVarListHandler *xVars, int *xMap, int *xPos);

	void iterateYVariables(TVarListHandler *xVars, int xId, int *yPos, int *yMin, int *yMax, int d);



};



class TShieldGeneratorTreeBase : public TShieldGeneratorBase {
public:

	// basic info
	THierarchicalPartition *yPartition;
	int dim; // dimension of space
	int lBottom,lTop; // specify which are the finest and coarsest relevant layer in yPartition
	double **yPos, **yRadii; // position and radii of y partition nodes
	int nLayers;
	
	// hierarchical data
	double **xposH;
	int *xresH;
	TVarListHandler **xNeighboursH;

	// data at current layer
	double *xpos; // position of x nodes
	int xres; // number of x nodes
	TVarListHandler *xNeighbours;


	TShieldGeneratorTreeBase(int _dim, THierarchicalPartition *_yPartition, double **_yPos, double **_yRadii,
			int _lBottom, int _lTop,
			int *_xresH, double **_xposH, TVarListHandler **_xNeighboursH);

	virtual ~TShieldGeneratorTreeBase();

	int setLayer(int newLayer);


	virtual void generateShield(TVarListHandler *xVars, TVarListHandler *xSupport);

	void addVariables_Shields(TVarListHandler *xVars, int *xMap, int x);
	void addVariables_Polytopes(TVarListHandler *xVars, int *xMap, int x);

	void iterateYVariables(TVarListHandler *xVars, int *xMap, int xA, int lParent, int yBParent);

	virtual bool checkConditionPlane(int xA, int x, int l, int yB, int y);
	virtual bool checkCondition(int xA, int l, int yB, int *xMap);


};

class TShieldGeneratorTreeBase_Benchmark : public TShieldGeneratorTreeBase {
public:

	TShieldGeneratorTreeBase_Benchmark(int _dim, THierarchicalPartition *_yPartition, double **_yPos, double **_yRadii,
			int _lBottom, int _lTop,
			int *_xresH, double **_xposH, TVarListHandler **_xNeighboursH);


	virtual ~TShieldGeneratorTreeBase_Benchmark();

	virtual void generateShield(TVarListHandler *xVars, TVarListHandler *xSupport);
	virtual bool checkCondition(int xA, int l, int yB, int *xMap);

};



template<class TReferenceTree>
class TShieldGeneratorTree_SqrEuclideanPrototype : public TReferenceTree {
public:
	using TReferenceTree::dim;
	using TReferenceTree::xpos;
	using TReferenceTree::yPos;
	using TReferenceTree::yRadii;
	using TReferenceTree::lTop;
	using TReferenceTree::lBottom;
	using TReferenceTree::shieldingTolerance;

	TShieldGeneratorTree_SqrEuclideanPrototype(int _dim, THierarchicalPartition *_yPartition, double **_yPos, double **_yRadii,
			int _lBottom, int _lTop,
			int *_xresH, double **_xposH, TVarListHandler **_xNeighboursH);

	bool checkConditionPlane(int xA, int x, int l, int yB, int y);
};

typedef TShieldGeneratorTree_SqrEuclideanPrototype<TShieldGeneratorTreeBase> TShieldGeneratorTree_SqrEuclidean;

template<class TReferenceTree>
class TShieldGeneratorTree_SqrEuclideanNoisePrototype : public TReferenceTree {
public:
	using TReferenceTree::dim;
	using TReferenceTree::xpos;
	using TReferenceTree::yPos;
	using TReferenceTree::yRadii;
	using TReferenceTree::lTop;
	using TReferenceTree::lBottom;
	using TReferenceTree::shieldingTolerance;

	double **c;
	double eta, lambda;
	int yres;

	TShieldGeneratorTree_SqrEuclideanNoisePrototype(int _dim,
			THierarchicalPartition *_yPartition, double **_yPos, double **_yRadii,
			int _lBottom, int _lTop,
			int *_xresH, double **_xposH, TVarListHandler **_xNeighboursH,
			double **_c, double _eta, double _lambda);
	//~TShieldGeneratorSqrEuclideanTree();
	bool checkConditionPlane(int xA, int x, int l, int yB, int y);
};


typedef TShieldGeneratorTree_SqrEuclideanNoisePrototype<TShieldGeneratorTreeBase> TShieldGeneratorTree_SqrEuclideanNoise;


template<class TReferenceTree>
class TShieldGeneratorTree_PEuclideanPrototype : public TReferenceTree {
public:
	using TReferenceTree::dim;
	using TReferenceTree::xpos;
	using TReferenceTree::yPos;
	using TReferenceTree::yRadii;
	using TReferenceTree::lTop;
	using TReferenceTree::lBottom;
	using TReferenceTree::shieldingTolerance;

	double p;
	double slack;
	double *xAyB,*xAys,*xsyB,*xsys,*xAxs; // memory for vectors to store: v=y-xA, vB=yB-xA, a=x-xA
	TShieldGeneratorTree_PEuclideanPrototype(int _dim, THierarchicalPartition *_yPartition, double **_yPos, double **_yRadii,
			int _lBottom, int _lTop,
			int *_xresH, double **_xposH, TVarListHandler **_xNeighboursH,
			double _p, double _slack);



	~TShieldGeneratorTree_PEuclideanPrototype();
	bool checkConditionPlane(int xA, int x, int l, int yB, int y);
	inline double getSubgradxAxs(double *z) {
		double innerProd,norm;
		norm=EUCL_innerProduct(z,z,dim);
		if (norm==0) {
			return 0;
		}

		innerProd=EUCL_innerProduct(z,xAxs,dim);

		return p*pow(norm,(p/2.-1))*innerProd;

	}
	inline double getPhi(double *z) {
		return pow(EUCL_innerProduct(z,z,dim),p/2);
	}
};

typedef TShieldGeneratorTree_PEuclideanPrototype<TShieldGeneratorTreeBase> TShieldGeneratorTree_PEuclidean;







#endif
