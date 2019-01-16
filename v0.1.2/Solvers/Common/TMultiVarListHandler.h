#ifndef TMULTIVARLISTHANDLER_H_
#define TMULTIVARLISTHANDLER_H_

#include<stdlib.h>
#include<vector>
#include<cmath>
#include"THierarchicalPartition.h"
#include<iostream>

using namespace std;

const double COST_INFINITY=1E10;

inline double getWFLogCost(double dist, double delta, double prefactor) {
	// post-processing for OT-FR cost
	
	if(dist>=M_PI*delta) {
		return COST_INFINITY;
	} else {
		return -2*log(cos(dist/(2*delta)))*prefactor;
	}
}

template <class T>
class TMultiVarListHandler {
public:
	int res,total,dim;
	vector<int> *lenList;
	vector<int*> **varList;
	vector<T> **signalList;
	TMultiVarListHandler(int _dim);
	TMultiVarListHandler(int _dim, int _res);
	virtual ~TMultiVarListHandler();
	void clear();
	void setupEmpty(int _res);
	void fillFromCSRIndexList(T *signal, int *indices, int *indptr, int _res, int _total);
	void writeToCSRIndexList(T *signal, int *indices, int *indptr);
	void addToLine(int x, T signal, int *yCandidate);
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class TMultiCostFunctionProvider {
public:
	double ***pos;
	double ***radii;
	double ***alpha;
	int dim, posDim, layerBottom;
	TMultiCostFunctionProvider(double ***_pos, double ***_radii, int _dim, int _posDim, int _layerBottom, double *** _alpha);
	virtual ~TMultiCostFunctionProvider();
	
	virtual inline double getCost(__attribute__((unused)) int layer, __attribute__((unused)) int *x) {
		return 0;
	}

	inline double getCostEff(int layer, int *x) {
		double costClean;
		int d;
		costClean=getCost(layer,x);
		// if no duals are given, just return standard cost
		if (alpha==NULL) { return costClean; }
		// go through all dual potentials and subtract them
		for(d=0;d<dim;d++) {
			costClean-=alpha[d][layer][x[d]];
		}
	
		return costClean;
	}
	
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

	static inline double EUCL_lincombSqr(double *a, double *b, double sa, double sb, int n) {
		// returns |sa*a+sb*b|^2, dimension given by n
		int i;
		double result,buffer;
		result=0;
		for(i=0;i<n;i++) {
			buffer=sa*a[i]+sb*b[i];
			result+=buffer*buffer;
		}
		return result;
	}

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Squared Euclidean */

class TMultiCostFunctionProvider_SquaredEuclidean : public TMultiCostFunctionProvider {
public:

	double weight;
	
	TMultiCostFunctionProvider_SquaredEuclidean(double ***_pos, double ***_radii, int _dim, int _posDim, int _layerBottom,
			double *** _alpha, double _weight);
			
	// only for dim=2. returns standard squared euclidean distance
	double getCost(int layer, int *x);
};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Squared Euclidean + Fisher-Rao */
class TMultiCostFunctionProvider_SquaredEuclideanWF : public TMultiCostFunctionProvider {
public:
	double delta; // weight for growth cost
	double cMax; // cuf-off for cost function value when approaching cut-locus
	double prefactor;

	TMultiCostFunctionProvider_SquaredEuclideanWF(double ***_pos, double ***_radii, int _dim, int _posDim, int _layerBottom,
			double *** _alpha, double _delta, double _cMax);
			
	// only for dim=2. returns standard squared euclidean distance
	double getCost(int layer, int *x);
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Squared Euclidean Barycenter */
class TMultiCostFunctionProvider_SquaredEuclideanBarycenter : public TMultiCostFunctionProvider {
public:
	double *lambda;
	TMultiCostFunctionProvider_SquaredEuclideanBarycenter(double ***_pos, double ***_radii, double *_lambda,
			int _dim, int _posDim, int _layerBottom, double *** _alpha);
	
	double getCost(int layer, int *x);

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Coulomb Cost */
class TMultiCostFunctionProvider_Coulomb : public TMultiCostFunctionProvider {
public:
	double *charges;
	TMultiCostFunctionProvider_Coulomb(double ***_pos, double ***_radii,
			int _dim, int _posDim, int _layerBottom, double *** _alpha,
			double *_charges);
	
	double getCost(int layer, int *x);

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Spherical Reflector */


class TMultiCostFunctionProvider_Reflector_Spherical : public TMultiCostFunctionProvider {
public:

	TMultiCostFunctionProvider_Reflector_Spherical(double ***_pos, double ***_radii, int _dim, int _posDim, int _layerBottom,
			double *** _alpha);
			
	// only for dim=2. returns spherical reflector cost
	double getCost(int layer, int *x);
};


const int CFPT_SquaredEuclidean=1;
const int CFPT_SquaredEuclideanBarycenter=2;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* auxiliary cost function providers for smoother level refinement */

/* interpolates cost of a cell with parent cell cost. used for smooth refinement */
class TMultiCostFunctionProvider_Interpolator : public TMultiCostFunctionProvider {
public:
	THierarchicalPartition **partition;
	TMultiCostFunctionProvider *coarse, *fine;
	double q;
	bool destroyChildren;
	static const double DBL_INFINITY=1E100;

	TMultiCostFunctionProvider_Interpolator(TMultiCostFunctionProvider *_coarse, TMultiCostFunctionProvider *_fine,
			THierarchicalPartition **_partition, double _q, double *** _alpha, bool _destroyChildren);
		
	~TMultiCostFunctionProvider_Interpolator();
	
	double getCost(int layer, int *x);

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* color transport */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Squared Euclidean + RGB cyclic */

class TMultiCostFunctionProvider_Color_SquaredEuclidean_RGB : public TMultiCostFunctionProvider {
public:

	double colorWeight;
	double FR_delta, FR_prefactor; // factors for Fisher-Rao mode
	bool FR_mode;
	int lTop;
	TMultiCostFunctionProvider_Color_SquaredEuclidean_RGB(double ***_pos, double ***_radii, double *** _alpha,
			int _dim, int _posDim, int _layerBottom,
			int _lTop, double _colorWeight, bool _FR_mode, double _FR_delta);
			
	// only for dim=2. returns standard squared euclidean distance
	double getCost(int layer, int *x);
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Squared Euclidean + HSV_HS */

class TMultiCostFunctionProvider_Color_SquaredEuclidean_HSV_HS : public TMultiCostFunctionProvider {
public:

	int liftMode;
	double colorWeight;
	double FR_delta, FR_prefactor; // factors for Fisher-Rao mode
	bool FR_mode;
	int lTop;
	double ***radiiHue, ***radiiVal;
	TMultiCostFunctionProvider_Color_SquaredEuclidean_HSV_HS(double ***_pos,
			double ***_radiiSpace, double ***_radiiHue, double ***_radiiVal,
			double *** _alpha,
			int _dim, int _posDim, int _layerBottom,
			int _liftMode,double _colorWeight, bool _FR_mode, double _FR_delta);
			
	// only for dim=2. returns standard squared euclidean distance
	double getCost(int layer, int *x);
};

#endif
