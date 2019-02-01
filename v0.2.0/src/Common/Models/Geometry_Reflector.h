#ifndef Model_Geometry_Reflector_H_
#define Model_Geometry_Reflector_H_

#include<Common/ErrorCodes.h>
#include<Common/Tools.h>
#include<Common/MultiScaleTools.h>
#include<Common/TCostFunctionProvider.h>
#include<Common/THierarchicalCostFunctionProvider.h>
#include<Common/Verbose.h>

// contains code for solving spherical reflector transport problems

// TCostFunctionProvider
// THierarchicalCostFunctionProvider

///////////////////////////////////////////////////////////////////////////////////////////////////	
// CostFunctionProvider
///////////////////////////////////////////////////////////////////////////////////////////////////	

class TCostFunctionProvider_Reflector : public TCostFunctionProvider_Dynamic {
public:

	TCostFunctionProvider_Reflector(int *_xresH, int *_yresH,
			double **_posXH, double **_posYH,
			int _nLayers, int _dim) :
					TCostFunctionProvider_Dynamic(_xresH,_yresH,_posXH,_posYH,_nLayers,_dim) {
	}		


	virtual inline double getCValue(int x, int y) {		
		double result;
		result=EUCL_innerProduct(posX+x*dim, posY+y*dim, dim);
		if(result>=1.) { return DBL_INFINITY; }
		result=-std::log(1-result);
		return result;
	}
};



///////////////////////////////////////////////////////////////////////////////////////////////////	
// HierarchicalCostFunctionProvider
///////////////////////////////////////////////////////////////////////////////////////////////////	


class THierarchicalCostFunctionProvider_Reflector : public THierarchicalCostFunctionProvider {
public:

	THierarchicalCostFunctionProvider_Reflector(
			double **_xPos, double **_yPos,
			double **_xRadii, double **_yRadii,
			int _posDim, int _layerBottom,
			bool _haveDuals,
			double **_alpha, double **_beta
			) :
			THierarchicalCostFunctionProvider(
					_xPos, _yPos,
					_xRadii, _yRadii,
					_posDim, _layerBottom,
					_haveDuals,
					_alpha, _beta) {
	};

	double getCostAsym(int layerX, int x, int layerY, int y);
};


#endif



