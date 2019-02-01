#ifndef TShieldGenerator_Models_H
#define TShieldGenerator_Models_H


#include<Common/ErrorCodes.h>
#include<Common/Tools.h>
#include<Common/Verbose.h>
#include<ShortCutSolver/TShieldGenerator.h>


/////////////////////////////////////////////////////////////////////////////////////////////
// metric^p on sphere
/////////////////////////////////////////////////////////////////////////////////////////////

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
			int *_xresH, double **_xposH, TVarListHandler **_xNeighboursH,
			double _p);


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

/////////////////////////////////////////////////////////////////////////////////////////////
// spherical reflector problem
/////////////////////////////////////////////////////////////////////////////////////////////


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
			int *_xresH, double **_xposH, TVarListHandler **_xNeighboursH);


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


#endif

