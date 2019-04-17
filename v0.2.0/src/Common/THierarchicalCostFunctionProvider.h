#ifndef THierarchicalCostFunctionProvider_H_
#define THierarchicalCostFunctionProvider_H_

#include<cstdlib>
#include<cmath>
#include<algorithm>

#include<Common/Tools.h>

class THierarchicalCostFunctionProvider {
public:
	static constexpr double DBL_INFINITY=1E100; // effective value for infinity

	double **xPos, **yPos; // pointers to coordinates of points:
	// hierarchical: xPos is list of pointers to coordinates at each hierarchy level
	double **xRadii, **yRadii; // likewise: radii of each hierarchical cell, used to compute lower bounds
	double **alpha, **beta; // pointers to hierarchical dual variables, to compute effective costs if required
	bool haveDuals; // indicates whether hierarchical dual variables are available
	int posDim; // dimensionality of coordinates in xPos and yPos arrays
	int layerBottom; // number of finest layer (0 is coarsest)
	
	THierarchicalCostFunctionProvider(
		double **_xPos, double **_yPos,
		double **_xRadii, double **_yRadii,
		int _posDim, int _layerBottom,
		bool _haveDuals,
		double **_alpha, double **_beta);
	
	virtual ~THierarchicalCostFunctionProvider();

	virtual void setLayerBottom(int _layerBottom);	

	virtual double getCost(int layer, int x, int y) {
		return getCostAsym(layer,x,layer,y);
	}

	virtual inline double getCostAsym(
			__attribute__((unused)) int layerX, __attribute__((unused)) int x,
			__attribute__((unused)) int layerY, __attribute__((unused)) int y) {
		return 0;
	}

	inline double getCostEff(int layer, int x, int y) {
		// compute effective cost, where dual variable values are subtracted
		if(haveDuals) {
			return getCost(layer,x,y)-alpha[layer][x]-beta[layer][y];
		} else {
			// if no duals available, return default cost
			return getCost(layer,x,y);
		}
	}

	inline double getCostEffAsym(int layerX, int x, int layerY, int y) {
		// compute effective cost, where dual variable values are subtracted
		if(haveDuals) {
			return getCostAsym(layerX,x,layerY,y)-alpha[layerX][x]-beta[layerY][y];
		} else {
			// if no duals available, return default cost
			return getCostAsym(layerX,x,layerY,y);
		}
	}

};

/* squared Euclidean distance in \R^n */

class THierarchicalCostFunctionProvider_SquaredEuclidean : public THierarchicalCostFunctionProvider {
public:

	double weight; // global rescaling parameter for Euclidean distance
	bool WFmode; // whether to compute cost function for Wasserstein--Fisher--Rao distance instead
	double WFlenscale; // max transport distance in WF mode
	double WFprefactor; // weight prefactor to avoid repeated computation
	
	THierarchicalCostFunctionProvider_SquaredEuclidean(
		double **_xPos, double **_yPos,
		double **_xRadii, double **_yRadii,
		int _posDim, int _layerBottom,
		bool _haveDuals,
		double **_alpha, double **_beta,
		double _weight,
		bool _WFmode,
		double _WFlenscale);

	THierarchicalCostFunctionProvider_SquaredEuclidean(
		double **_xPos, double **_yPos,
		double **_xRadii, double **_yRadii,
		int _posDim, int _layerBottom,
		bool _haveDuals,
		double **_alpha, double **_beta,
		double _weight) : 

		THierarchicalCostFunctionProvider_SquaredEuclidean(
			_xPos, _yPos,
			_xRadii, _yRadii,
			_posDim, _layerBottom,
			_haveDuals,
			_alpha, _beta,
			_weight,
			false,
			0) {};

	void setWFlenscale(const double _WFlenscale);
	double getCostAsym(int layerX, int x, int layerY, int y);
};


/* Wp: Euclidean distance in \R^n to power p*/

class THierarchicalCostFunctionProvider_PEuclidean : public THierarchicalCostFunctionProvider {
public:

	double weight; // global rescaling parameter for Euclidean distance
	double p; // exponent for distance
	
	THierarchicalCostFunctionProvider_PEuclidean(
		double **_xPos, double **_yPos,
		double **_xRadii, double **_yRadii,
		int _posDim, int _layerBottom,
		bool _haveDuals,
		double **_alpha, double **_beta,
		double _weight,
		double _p
		);

	double getCostAsym(int layerX, int x, int layerY, int y);
};


/* Wp: Euclidean distance in \R^n to power p*/

class THierarchicalCostFunctionProvider_Hyperbolic : public THierarchicalCostFunctionProvider {
public:

	double scale; // scale for hyperboloid: x_0^2 = scale^2 + \sum_{i=1}^n x_i^2
	double scaleSqr; // to avoid computing the square all the time
	
	// in this implementation, pos only stores coordinates x_1 to x_n; x_0 is implicit via above "Minkowski" formula
	
	THierarchicalCostFunctionProvider_Hyperbolic(
		double **_xPos, double **_yPos,
		double **_xRadii, double **_yRadii,
		int _posDim, int _layerBottom,
		bool _haveDuals,
		double **_alpha, double **_beta,
		double _scale
		);

	double getCostAsym(int layerX, int x, int layerY, int y);
};


#endif
