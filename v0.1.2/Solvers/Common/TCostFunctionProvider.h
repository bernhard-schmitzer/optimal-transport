#ifndef TCostFunctionProvider_H_
#define TCostFunctionProvider_H_

#include<stdlib.h>
#include<vector>
#include"TVarListHandler.h"

using namespace std;

class TCostFunctionProviderBase {
	public:
	
	virtual ~TCostFunctionProviderBase() {
	}

	virtual bool free_c() {
		return false;
	}
	virtual double* getC(__attribute__((unused)) TVarListHandler *xVars) {
		return NULL;
	}
};

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
			for(yIndex=0;yIndex<xVars->lenList->at(x);yIndex++) {
				y=xVars->varList[x]->at(yIndex);
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

#endif
