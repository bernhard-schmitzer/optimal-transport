#ifndef TShieldGenerator_H
#define TShieldGenerator_H

#include<TVarListHandler.h>
#include<THierarchicalPartition.h>
#include<cmath>
#include<iostream>

using namespace std;


void GridToolsGetStrides(int dim, int *dims, int *strides);
int GridToolsGetIdFromPos(int dim, int *pos, int *strides);
void GridToolsGetPosFromId(int dim, int id, int *pos, int *strides);

void GridToolsGetNeighbours(int dim, int *dims, TVarListHandler *neighbours);
void GridToolsGetNeighbours_Torus(int dim, int *dims, int torusDim, TVarListHandler *neighbours);
void GridToolsGetNeighbours_Torus_iterateXVariables(TVarListHandler *neighbours, int *pos, int *dims, int *strides,
		int dim, int torusDim, int d);



class TShieldGeneratorBase {
public:
	static const double shieldingTolerance=1E-5;
	// benchmark variable to be used by some shielding generators
	int n_shielding_queries;
	TShieldGeneratorBase() { n_shielding_queries=0; };
	virtual ~TShieldGeneratorBase() {};
	virtual void generateShield(__attribute__((unused)) TVarListHandler *xVars,
			__attribute__((unused)) TVarListHandler *xSupport) {};
	static void getXMap(int *xMap, TVarListHandler *xSupport);
};


class TShieldGeneratorGrid_SqrEuclidean : public TShieldGeneratorBase {
public:
	int *xDims,*yDims,dim; // dim: dimension of space, xDims,yDims: number of points along each axis
	int *xStrides, *yStrides; // strides: difference in index, for moving on grid point in a given dimension
	TShieldGeneratorGrid_SqrEuclidean(int _dim, int *_xDims, int *_yDims);
	~TShieldGeneratorGrid_SqrEuclidean();

	void generateShield(TVarListHandler *xVars, TVarListHandler *xSupport);
	void iterateXVariables(TVarListHandler *xVars, int *xMap, int *xPos, int d);

	void addVariables_Shields(TVarListHandler *xVars, int *xMap, int *xPos);
	void addVariables_Rectangles(TVarListHandler *xVars, int *xMap, int *xPos);

	void iterateYVariables(TVarListHandler *xVars, int xId, int *yPos, int *yMin, int *yMax, int d);

};

/*
class TShieldGeneratorGrid_Torus : public TShieldGeneratorBase {
public:
	int *xDims,*yDims,dim; // dim: dimension of space, xDims,yDims: number of points along each axis
	int *xStrides, *yStrides; // strides: difference in index, for moving on grid point in a given dimension
	double *radii;
	TShieldGeneratorGrid_Torus(int _dim, int *_xDims, int *_yDims, double *_radii);
	~TShieldGeneratorGrid_Torus();

	void generateShield(TVarListHandler *xVars, TVarListHandler *xSupport);
	void iterateXVariables(TVarListHandler *xVars, int *xMap, int *xPos, int d);

	void addVariables_Shields(TVarListHandler *xVars, int *xMap, int *xPos);
	void addVariables_Rectangles(TVarListHandler *xVars, int *xMap, int *xPos);

	void iterateYVariables(TVarListHandler *xVars, int xId, int *yPos, int *yMin, int *yMax, int d);

};
*/


class TShieldGeneratorGrid_Padding : public TShieldGeneratorBase {
public:
	int *xDims,*yDims,dim; // dim: dimension of space, xDims,yDims: number of points along each axis
	int *xStrides, *yStrides; // strides: difference in index, for moving on grid point in a given dimension
	int width; // width of padding hull. for now: simply call inner padding routine width times.
	TShieldGeneratorGrid_Padding(int _dim, int *_xDims, int *_yDims, int _width);
	~TShieldGeneratorGrid_Padding();

	void generateShield(TVarListHandler *xVars, TVarListHandler *xSupport);
	void iterateVariables(TVarListHandler *xVars, TVarListHandler *xSupport);

	void addVariables(TVarListHandler *xVars, int xId, int yId, int *xPos, int *yPos);

};


class TShieldGeneratorTreeBase : public TShieldGeneratorBase {
public:
	int dim; // dimension of space
	THierarchicalPartition *yPartition;
	double **yPos, **yRadii; // position and radii of y partition nodes
	int lBottom,lTop; // specify which are the finest and coarsest relevant layer in yPartition
	double *xpos; // position of x nodes
	int xres; // number of x nodes
	TVarListHandler *xNeighbours;

	// Benchmark Variables
	//int *nCalls_checkCondition, *nCalls_checkConditionPlane;


	TShieldGeneratorTreeBase(int _dim, THierarchicalPartition *_yPartition, double **_yPos, double **_yRadii,
			int _lBottom, int _lTop,
			double *_xpos, TVarListHandler *_xNeighbours);

	virtual ~TShieldGeneratorTreeBase();

	virtual void generateShield(TVarListHandler *xVars, TVarListHandler *xSupport);

	void addVariables_Shields(TVarListHandler *xVars, int *xMap, int x);
	void addVariables_Polytopes(TVarListHandler *xVars, int *xMap, int x);

	void iterateYVariables(TVarListHandler *xVars, int *xMap, int xA, int lParent, int yBParent);

	virtual bool checkConditionPlane(int xA, int x, int l, int yB, int y);
	virtual bool checkCondition(int xA, int l, int yB, int *xMap);

	static inline double EUCL_innerProduct(double *a, double *b, int n) {
		// computes inner product a.b, dimension given by n
		double result=0;
		int i;
		for(i=0;i<n;i++) {
			result+=a[i]*b[i];
		}
		return result;
	}

	static inline void EUCL_lincomb(double *a, double *b, double *c, double sa, double sb, int n) {
		// stores sa*a+sb*b to c, dimension given by n
		int i;
		for(i=0;i<n;i++) {
			c[i]=sa*a[i]+sb*b[i];
		}
	}


};


class TShieldGeneratorTreeBase_Benchmark : public TShieldGeneratorTreeBase {
public:

	TShieldGeneratorTreeBase_Benchmark(int _dim, THierarchicalPartition *_yPartition, double **_yPos, double **_yRadii,
			int _lBottom, int _lTop,
			double *_xpos, TVarListHandler *_xNeighbours);

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
			double *_xpos, TVarListHandler *_xNeighbours);
	//~TShieldGeneratorSqrEuclideanTree();
	bool checkConditionPlane(int xA, int x, int l, int yB, int y);
};

typedef TShieldGeneratorTree_SqrEuclideanPrototype<TShieldGeneratorTreeBase> TShieldGeneratorTree_SqrEuclidean;


template<class TReferenceTree>
class TShieldGeneratorTree_TorusPrototype : public TReferenceTree {
public:
	using TReferenceTree::dim;
	using TReferenceTree::xpos;
	using TReferenceTree::yPos;
	using TReferenceTree::yRadii;
	using TReferenceTree::lTop;
	using TReferenceTree::lBottom;
	using TReferenceTree::shieldingTolerance;
	using TReferenceTree::xNeighbours;

	double ***yTorusRadii;
	double *torusRadii;
	int torusDim;

	TShieldGeneratorTree_TorusPrototype(int _dim, THierarchicalPartition *_yPartition, double **_yPos, double **_yRadii,
			double ***_yTorusRadii,
			int _lBottom, int _lTop,
			double *_xpos, TVarListHandler *_xNeighbours,
			double *_torusRadii, int _torusDim);
	bool checkCondition(int xA, int l, int yB, int *xMap);
	double slackConditionPlane(int xA, int x, int l, int yB, int y);
	double slackConditionS1(int xA, int x, int l, int yB, int y, int axis);

	inline double translateCoord(double coord, double r, double delta) {
		double result;
		result=(coord/r+delta);
		while(result<0) {
			result+=1;
		}
		while(result>1) {
			result-=1;
		}
		return result;
	}

	inline double psi(double y, double xs, double rad) {
		if(y<xs-0.5) {
			return min(y+rad,xs-0.5)*(2*xs-3)+0.25-pow(xs-1,2);
		} else {
			return max(y-rad,xs-0.5)*(2*xs-1)+0.25-pow(xs,2);
		}
	}
};

typedef TShieldGeneratorTree_TorusPrototype<TShieldGeneratorTreeBase> TShieldGeneratorTree_Torus;


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
			double *_xpos, TVarListHandler *_xNeighbours,
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
			double *_xpos, TVarListHandler *_xNeighbours, double _p, double _slack);
	~TShieldGeneratorTree_PEuclideanPrototype();
	bool checkConditionPlane(int xA, int x, int l, int yB, int y);
	inline double getSubgradxAxs(double *z) {
		double innerProd,norm;
		norm=TReferenceTree::EUCL_innerProduct(z,z,dim);
		if (norm==0) {
			return 0;
		}

		innerProd=TReferenceTree::EUCL_innerProduct(z,xAxs,dim);

		return p*pow(norm,(p/2.-1))*innerProd;

	}
	inline double getPhi(double *z) {
		return pow(TReferenceTree::EUCL_innerProduct(z,z,dim),p/2);
	}
};

typedef TShieldGeneratorTree_PEuclideanPrototype<TShieldGeneratorTreeBase> TShieldGeneratorTree_PEuclidean;


template<class TReferenceTree>
class TShieldGeneratorTree_SpherePrototype : public TReferenceTree {
public:
	using TReferenceTree::dim;
	using TReferenceTree::xpos;
	using TReferenceTree::yPos;
	using TReferenceTree::yRadii;
	using TReferenceTree::lTop;
	using TReferenceTree::lBottom;
	using TReferenceTree::shieldingTolerance;
	double p;
	TShieldGeneratorTree_SpherePrototype(int _dim, THierarchicalPartition *_yPartition, double **_yPos, double **_yRadii,
			int _lBottom, int _lTop,
			double *_xpos, TVarListHandler *_xNeighbours, double _p);
	~TShieldGeneratorTree_SpherePrototype();
	bool checkConditionPlane(int xA, int x, int l, int yB, int y);

	inline double getSubgrad(double z) {
		return p*pow(z,(p-1));
	}

	inline double getPhi(double z) {
		return pow(z,p);
	}
};

typedef TShieldGeneratorTree_SpherePrototype<TShieldGeneratorTreeBase> TShieldGeneratorTree_Sphere;


template<class TReferenceTree>
class TShieldGeneratorTree_ReflectorPrototype : public TReferenceTree {
public:
	using TReferenceTree::dim;
	using TReferenceTree::xpos;
	using TReferenceTree::yPos;
	using TReferenceTree::yRadii;
	using TReferenceTree::lTop;
	using TReferenceTree::lBottom;
	using TReferenceTree::shieldingTolerance;
	TShieldGeneratorTree_ReflectorPrototype(int _dim, THierarchicalPartition *_yPartition, double **_yPos, double **_yRadii,
			int _lBottom, int _lTop,
			double *_xpos, TVarListHandler *_xNeighbours);
	~TShieldGeneratorTree_ReflectorPrototype();
	bool checkConditionPlane(int xA, int x, int l, int yB, int y);

	inline double getSubgrad(double z) {
		return -sin(z)/(1-cos(z));
	}

	inline double getPhi(double z) {
		return -log(1-cos(z));
	}
};

typedef TShieldGeneratorTree_ReflectorPrototype<TShieldGeneratorTreeBase> TShieldGeneratorTree_Reflector;

class TShieldGeneratorComparison : public TShieldGeneratorBase {
public:
	TShieldGeneratorBase *shieldGenerator1, *shieldGenerator2;
	TShieldGeneratorComparison(TShieldGeneratorBase *_shieldGenerator1, TShieldGeneratorBase *_shieldGenerator2);
	void generateShield(TVarListHandler *xVars, TVarListHandler *xSupport);
};



class TShieldingVerification {
public:
	double *c;
	int xres,yres;
	TVarListHandler *xNeighbours;
	TShieldingVerification(double *_c, int _xres, int _yres, TVarListHandler *_xNeighbours);
	TVarListHandler* verify(TVarListHandler *xVars, int *xMap);
};


class TShieldingVerificationDuplex {
public:
	double *c;
	int xres,yres;
	TVarListHandler *xNeighbours;
	TVarListHandler *yNeighbours;
	TShieldingVerificationDuplex(double *_c, int _xres, int _yres,
			TVarListHandler *_xNeighbours, TVarListHandler *_yNeighbours);
	TVarListHandler* verify(TVarListHandler *xVars, TVarListHandler *xSupport, TVarListHandler *ySupport);
};


#endif
