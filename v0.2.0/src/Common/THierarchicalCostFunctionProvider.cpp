#include"THierarchicalCostFunctionProvider.h"

THierarchicalCostFunctionProvider::THierarchicalCostFunctionProvider(
		double **_xPos, double **_yPos,
		double **_xRadii, double **_yRadii,
		int _posDim, int _layerBottom,
		bool _haveDuals,
		double **_alpha, double **_beta) {
	xPos=_xPos;
	yPos=_yPos;
	xRadii=_xRadii;
	yRadii=_yRadii;
	posDim=_posDim;
	layerBottom=_layerBottom;
	haveDuals=_haveDuals;
	if(haveDuals) {
		alpha=_alpha;
		beta=_beta;
	} else {
		alpha=NULL;
		beta=NULL;
	}
}


THierarchicalCostFunctionProvider::~THierarchicalCostFunctionProvider() {
}


void THierarchicalCostFunctionProvider::setLayerBottom(int _layerBottom) {
	layerBottom=_layerBottom;	
}


/* squared Euclidean distance */

THierarchicalCostFunctionProvider_SquaredEuclidean::THierarchicalCostFunctionProvider_SquaredEuclidean(
		double **_xPos, double **_yPos,
		double **_xRadii, double **_yRadii,
		int _posDim, int _layerBottom,
		bool _haveDuals,
		double **_alpha, double **_beta,
		double _weight,
		bool _WFmode,
		double _WFlenscale) :
			THierarchicalCostFunctionProvider(
				_xPos, _yPos,
				_xRadii, _yRadii,
				_posDim, _layerBottom,
				_haveDuals,
				_alpha, _beta) {
	weight=_weight;
	WFmode=_WFmode;
	if(WFmode) {
		WFlenscale=_WFlenscale;
		WFprefactor=4*WFlenscale*WFlenscale;
	} else {
		WFlenscale=0;
		WFprefactor=0;
	}
	
}

void THierarchicalCostFunctionProvider_SquaredEuclidean::setWFlenscale(const double _WFlenscale) {
	WFlenscale=_WFlenscale;
	WFprefactor=4*WFlenscale*WFlenscale;
}


double THierarchicalCostFunctionProvider_SquaredEuclidean::getCostAsym(int layerX, int x, int layerY, int y) {
	double result;
	
	// compute squared Euclidean distance between representatives
	result=EUCL_lincombSqr(xPos[layerX]+(x*posDim), yPos[layerY]+(y*posDim), 1, -1, posDim);
	result=std::sqrt(result);

	if((layerX<layerBottom) || (layerY<layerBottom)) {
		// if not at finest layer, need to compute lower bound
		if(layerX<layerBottom) {
			result-=xRadii[layerX][x];
		}
		if(layerY<layerBottom) {
			result-=yRadii[layerY][y];
		}
		if(result<0) {
			result=0;
		}
	}

	result=weight*result;

	if(!WFmode) {
		return std::pow(result,2);
	} else {
		if(result>M_PI*WFlenscale) {
			return DBL_INFINITY;
		} else {
			return WFprefactor*(-2*std::log(std::cos(result/(2*WFlenscale))));
		}
		
	}

}



/* Euclidean distance to power p */

THierarchicalCostFunctionProvider_PEuclidean::THierarchicalCostFunctionProvider_PEuclidean(
		double **_xPos, double **_yPos,
		double **_xRadii, double **_yRadii,
		int _posDim, int _layerBottom,
		bool _haveDuals,
		double **_alpha, double **_beta,
		double _weight,
		double _p) :
			THierarchicalCostFunctionProvider(
				_xPos, _yPos,
				_xRadii, _yRadii,
				_posDim, _layerBottom,
				_haveDuals,
				_alpha, _beta) {
	weight=_weight;
	p=_p;
	
}


double THierarchicalCostFunctionProvider_PEuclidean::getCostAsym(int layerX, int x, int layerY, int y) {
	double result;
	
	// compute squared Euclidean distance between representatives, fixed to dim=2
	result=EUCL_lincombSqr(xPos[layerX]+(x*posDim), yPos[layerY]+(y*posDim), 1, -1, posDim);
	result=sqrt(result);

	if((layerX<layerBottom) || (layerY<layerBottom)) {
		// if not at finest layer, need to compute lower bound
		if(layerX<layerBottom) {
			result-=xRadii[layerX][x];
		}
		if(layerY<layerBottom) {
			result-=yRadii[layerY][y];
		}
		if(result<0) {
			result=0;
		}
	}

	result=weight*result;
	result=pow(result,p);

	return result;
}


/* hyperbolic space */


THierarchicalCostFunctionProvider_Hyperbolic::THierarchicalCostFunctionProvider_Hyperbolic(
		double **_xPos, double **_yPos,
		double **_xRadii, double **_yRadii,
		int _posDim, int _layerBottom,
		bool _haveDuals,
		double **_alpha, double **_beta,
		double _scale
		) :
			THierarchicalCostFunctionProvider(
				_xPos, _yPos,
				_xRadii, _yRadii,
				_posDim, _layerBottom,
				_haveDuals,
				_alpha, _beta) {
	scale=_scale;
	scaleSqr=std::pow(scale,2);
	
}


double THierarchicalCostFunctionProvider_Hyperbolic::getCostAsym(int layerX, int x, int layerY, int y) {
	// (u,v): points on hyperboloid. uHat, vHat: only "spatial" coordinates 1..n (the ones, stores in pos)
	// distance formula: d(u,v) = arccosh(\sqrt{(1+uHat^2) (1+vHat^2)}- <uHat,vHat>)
	// call argument of arccosh omega
	
	// careful: have added parameter scale. effectively, uHat stores real uHat*scale.
	// so for large scale: zoom in near apex of hyperboloid, large coordinates are still close to apex.
	
	// first compute uHat^2, vHat^2 and <uHat,vHat>:
	
	
//	for(int i=0;i<posDim;i++) {
//		printf("%f ",xPos[layerX][x*posDim+i]);
//	}
//	printf("\n");

//	for(int i=0;i<posDim;i++) {
//		printf("%f ",yPos[layerY][y*posDim+i]);
//	}
//	printf("\n");
	
	double uHatLen=std::sqrt(EUCL_innerProduct(xPos[layerX]+(x*posDim),xPos[layerX]+(x*posDim),posDim));
	double vHatLen=std::sqrt(EUCL_innerProduct(yPos[layerY]+(y*posDim),yPos[layerY]+(y*posDim),posDim));
	double uHatvHat=EUCL_innerProduct(xPos[layerX]+(x*posDim),yPos[layerY]+(y*posDim),posDim);

	double uHatSqr, vHatSqr;

	//printf("uHatLen: %f\nvHatLen: %f\nuHatvHat: %f\n",uHatLen,vHatLen,uHatvHat);

	
	// if not at finest layer, need to compute lower bound
	if(layerX<layerBottom) {
		//printf("xrad: %f\n",xRadii[layerX][x]);
		uHatSqr=std::pow(std::max(0.,uHatLen-xRadii[layerX][x]),2);
		uHatvHat+=vHatLen*xRadii[layerX][x];
	} else {
		uHatSqr=std::pow(uHatLen,2);
	}
	
	if(layerY<layerBottom) {
		//printf("yrad: %f\n",yRadii[layerY][y]);
		vHatSqr=std::pow(std::max(0.,vHatLen-yRadii[layerY][y]),2);
		uHatvHat+=uHatLen*yRadii[layerY][y];
	} else {
		vHatSqr=std::pow(vHatLen,2);
	}

	if((layerX<layerBottom) && (layerY<layerBottom)) {
		uHatvHat+=xRadii[layerX][x]*yRadii[layerY][y];
	}
	
	double omega=(std::sqrt((scaleSqr+uHatSqr)*(scaleSqr+vHatSqr))-uHatvHat)/scaleSqr;
	
	//printf("uHatvHat: %f\nuHatSqr: %f\nvHatSqr: %f\nomega: %f\nresult: %f\n",uHatvHat,uHatSqr,vHatSqr,omega,std::acosh(omega));
	
	if(omega<=1) {
		return 0;
	}
	
	return scaleSqr*std::pow(std::acosh(omega),2);


}

