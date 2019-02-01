#ifndef TCostFunctionProvider_H_
#define TCostFunctionProvider_H_

#include<cstdlib>
#include<Common/ErrorCodes.h>
#include<Common/TVarListHandler.h>

class TCostFunctionProviderBase {
	public:
	static constexpr double DBL_INFINITY=1E100; // effective value for infinity
	
	virtual ~TCostFunctionProviderBase() {
	}

	virtual bool free_c() {
		return false;
	}
	virtual double* getC(__attribute__((unused)) TVarListHandler *xVars) {
		return NULL;
	}
	
	virtual int setLayer(__attribute__((unused)) int newLayer) {
		eprintf("ERROR: called setLayer() on TCostFunctionProviderBase\n");
		return ERR_BASE_NOTIMPLEMENTED;
	}
};

// just returns a previously allocated array.
// only used when solving dense problems. not relevant in sparse and multi-scale algorithms
class TCostFunctionProviderStatic : public TCostFunctionProviderBase {
	public:
	double *c;
	TCostFunctionProviderStatic(double *_c) {
		c=_c;
	}
	double* getC(__attribute__((unused)) TVarListHandler *xVars) {
		return c;
	}
};

// pseudo sparse variant
// full cost is computed and stored. sparse subsets are returned on demand
class TCostFunctionProviderFullArray : public TCostFunctionProviderBase {
	public:
	double *c;
	int xres, yres;
	TCostFunctionProviderFullArray(int _xres, int _yres, double *_c) {
		xres=_xres;
		yres=_yres;
		c=_c;
	}
	double* getC(TVarListHandler *xVars) {
		double *result;
		result=(double*) malloc(sizeof(double)*xVars->total);
		int x,yIndex,y,offset;
		offset=0;
		for(x=0;x<xVars->res;x++) {
			for(yIndex=0;yIndex<xVars->lenList[x];yIndex++) {
				y=xVars->varList[x][yIndex];
				result[offset]=c[x*yres+y];
				offset++;
			}
		}
			
		return result;
	}
	virtual bool free_c() {
		return true;
	}

};


// prototype of truly dynamic version:
// cost of relevant entries is computed on demand
// in prototype version: squared Euclidean distance
class TCostFunctionProvider_Dynamic : public TCostFunctionProviderBase {
	public:

	// basic data
	int dim; // spatial dimension of coordinates
	int nLayers; // number of layers
	int layer; // number of currently active layer
	
	// hierarchical data over all layers
	int *xresH, *yresH;
	double **posXH, **posYH;

	// data at current layer
	int xres, yres;
	double *posX, *posY;
	
	TCostFunctionProvider_Dynamic(int *_xresH, int *_yresH,
		double **_posXH, double **_posYH,
		int _nLayers, int _dim) {
		
		xresH=_xresH;
		yresH=_yresH;
		posXH=_posXH;
		posYH=_posYH;
		nLayers=_nLayers;
		dim=_dim;
		
		xres=0;
		yres=0;
		posX=NULL;
		posY=NULL;
		setLayer(0);
		
	}
	
	int setLayer(int newLayer) {
		if(newLayer>=nLayers) {
			return ERR_MULTISCALE_EXCEEDEDLEVELS;
		}
		layer=newLayer;
		xres=xresH[layer];
		yres=yresH[layer];
		posX=posXH[layer];
		posY=posYH[layer];
		return 0;
	}
		
	double* getC(TVarListHandler *xVars) {
		double *result;
		result=(double*) malloc(sizeof(double)*xVars->total);
		int x,yIndex,y,offset;
		offset=0;
		for(x=0;x<xVars->res;x++) {
			for(yIndex=0;yIndex<xVars->lenList[x];yIndex++) {
				y=xVars->varList[x][yIndex];
				result[offset]=getCValue(x,y);
				offset++;
			}
		}
			
		return result;
	}
	
	// returns full xres*yres cost function array
	// useful aux function and for debugging
	void getCDense(double *c) {
		int x,y;
	 	for(x=0;x<xres;x++) {
	 		for(y=0;y<yres;y++) {
	 			c[x*yres+y]=getCValue(x,y);	 			
	 		}
	 	}
 	}

	
	virtual inline double getCValue(int x, int y) {		
		double result,diff;
		result=0;
		for(int i=0;i<dim;i++) {
			diff=posX[x*dim+i]-posY[y*dim+i];
			result+=diff*diff;
		}
		return result;
	}
	
	virtual bool free_c() {
		return true;
	}

};

class TCostFunctionProvider_PEuclidean : public TCostFunctionProvider_Dynamic {
public:

	double p;

	TCostFunctionProvider_PEuclidean(int *_xresH, int *_yresH,
			double **_posXH, double **_posYH,
			int _nLayers, int _dim, double _p) :
					TCostFunctionProvider_Dynamic(_xresH,_yresH,_posXH,_posYH,_nLayers,_dim) {
		p=_p;
	}		


	virtual inline double getCValue(int x, int y) {		
		double result,diff;
		result=0;
		for(int i=0;i<dim;i++) {
			diff=posX[x*dim+i]-posY[y*dim+i];
			result+=diff*diff;
		}
		return std::pow(result,p/2.);
	}
};


#endif
