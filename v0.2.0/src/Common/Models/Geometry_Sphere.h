#ifndef Model_Geometry_Sphere_H_
#define Model_Geometry_Sphere_H_

#include<Common/ErrorCodes.h>
#include<Common/Tools.h>
#include<Common/MultiScaleTools.h>
#include<Common/TCostFunctionProvider.h>
#include<Common/THierarchicalCostFunctionProvider.h>
#include<Common/Verbose.h>

// contains code for solving transport problems on the sphere

// TCostFunctionProvider for sphere distance (used primarily in ShortCutSolver)
// THierarchicalCostFunctionProvider for sphere distance (used primarily in SinkhornSolver)
// TMultiScaleSetupBase for setup on sphere
	// this inludes: projecting hierarchical tree points back onto sphere
	// computing radii in sphere metric

///////////////////////////////////////////////////////////////////////////////////////////////////	
// CostFunctionProvider
///////////////////////////////////////////////////////////////////////////////////////////////////	

class TCostFunctionProvider_Sphere : public TCostFunctionProvider_Dynamic {
public:

	double p;

	TCostFunctionProvider_Sphere(int *_xresH, int *_yresH,
			double **_posXH, double **_posYH,
			int _nLayers, int _dim, double _p) :
					TCostFunctionProvider_Dynamic(_xresH,_yresH,_posXH,_posYH,_nLayers,_dim) {
		p=_p;
	}		


	virtual inline double getCValue(int x, int y) {		
		double result;
		result=EUCL_innerProduct(posX+x*dim, posY+y*dim, dim);
		if(result>=1.) { return 0.; }
		if(result<=-1.) { return std::pow(M_PI,p); }
		
		result=std::acos(result);
		return std::pow(result,p);
	}
};


///////////////////////////////////////////////////////////////////////////////////////////////////	
// MultiScaleSetup
///////////////////////////////////////////////////////////////////////////////////////////////////	


class TMultiScaleSetupSingleSphere : public TMultiScaleSetupSingleBase {
public:
	
	static constexpr double sphereCenterTolerance=1E-14;
	// points with length smaller than this will be projected to northpole

	// inherit constructor
	TMultiScaleSetupSingleSphere(TDoubleMatrix *_pos, double *_mu, int _depth) :
		TMultiScaleSetupSingleBase(_pos,_mu,_depth) {};
	
	int SetupProjectPoints();
	int SetupProjectPoints_Array(double *pos, int n, int dim);
	int SetupRadii();
	
	double SphereDistance(double *posx, double *posy);
};



///////////////////////////////////////////////////////////////////////////////////////////////////	
// HierarchicalCostFunctionProvider
///////////////////////////////////////////////////////////////////////////////////////////////////	


class THierarchicalCostFunctionProvider_Sphere : public THierarchicalCostFunctionProvider {
public:

	double p;

	THierarchicalCostFunctionProvider_Sphere(
			double **_xPos, double **_yPos,
			double **_xRadii, double **_yRadii,
			int _posDim, int _layerBottom,
			bool _haveDuals,
			double **_alpha, double **_beta,
			double _p
			) :
			THierarchicalCostFunctionProvider(
					_xPos, _yPos,
					_xRadii, _yRadii,
					_posDim, _layerBottom,
					_haveDuals,
					_alpha, _beta) {
		p=_p;
	};

	double getCostAsym(int layerX, int x, int layerY, int y);
};


#endif



