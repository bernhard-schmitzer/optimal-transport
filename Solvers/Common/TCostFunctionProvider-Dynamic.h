#ifndef TCostFunctionProvider_Dynamic_H_
#define TCostFunctionProvider_Dynamic_H_

#include<stdlib.h>
#include<vector>
#include"TVarListHandler.h"
#include"TCostFunctionProvider.h"

using namespace std;

class TCostFunctionProvider_Dynamic : public TCostFunctionProviderBase {
	public:
	double *c;
	int xres, yres, dim;
	double *xPos, *yPos;
	
	TCostFunctionProvider_Dynamic(int _xres, int _yres,
		double *_xPos, double *_yPos,
		int _dim) {
		
		xres=_xres;
		yres=_yres;
		xPos=_xPos;
		yPos=_yPos;
		dim=_dim;
		
	}
		
	double* getC(TVarListHandler *xVars) {
		double *result;
		result=(double*) malloc(sizeof(double)*xVars->total);
		int x,yIndex,y,offset;
		offset=0;
		for(x=0;x<xVars->res;x++) {
			for(yIndex=0;yIndex<xVars->lenList->at(x);yIndex++) {
				y=xVars->varList[x]->at(yIndex);
				result[offset]=getCValue(x,y);
				offset++;
			}
		}
			
		return result;
	}
	
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
			diff=xPos[x*dim+i]-yPos[y*dim+i];
			result+=diff*diff;
		}
		return result;
	}
	
	virtual bool free_c() {
		return true;
	}

};


class TCostFunctionProvider_Dynamic_Torus : public TCostFunctionProvider_Dynamic {
	public:
	double *radius;
	int torusDim;
	
	TCostFunctionProvider_Dynamic_Torus(int _xres, int _yres,
		double *_xPos, double *_yPos,
		int _dim,
		double *_radius, int _torusDim) : TCostFunctionProvider_Dynamic(_xres,_yres,_xPos,_yPos,_dim) {

		radius=_radius;
		torusDim=_torusDim;
		
	}
		
	virtual inline double getCValue(int x, int y) {
		double result,diff;
		result=0;
		for(int i=0;i<dim;i++) {
			diff=abs(xPos[x*dim+i]-yPos[y*dim+i]);
			if(i<torusDim) {
				diff=min(diff,radius[i]-diff);
			}
			result+=diff*diff;
		}
		return result;
	}
	
};


#endif

