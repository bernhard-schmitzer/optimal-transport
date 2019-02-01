#include"Geometry_Reflector.h"

	
///////////////////////////////////////////////////////////////////////////////////////////////////	
// HierarchicalCostFunctionProvider
///////////////////////////////////////////////////////////////////////////////////////////////////	


double THierarchicalCostFunctionProvider_Reflector::getCostAsym(int layerX, int x, int layerY, int y) {
	double result;
	
	// compute sphere distance from inner product
	result=EUCL_innerProduct(xPos[layerX]+(x*posDim), yPos[layerY]+(y*posDim), posDim);
	if(result>=1.) {
		result=0.;
	} else {
		if(result<=-1.) {
			result=M_PI;
		} else {
			result=std::acos(result);
		}
	}

	
	// if not at finest layer, need to compute lower bound
	// for reflector this means, get upper bound on distance
	if(layerX<layerBottom) {
		result+=xRadii[layerX][x];
	}
	if(layerY<layerBottom) {
		result+=yRadii[layerY][y];
	}
	// truncate at M_PI.
	if(result>M_PI) {
		result=M_PI;
	}

	return -std::log(1-std::cos(result));
}


