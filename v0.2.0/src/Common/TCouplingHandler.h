#ifndef TCouplingHandler_H_
#define TCouplingHandler_H_

#include<cstdlib>
#include<cmath>
#include<Common/TVarListHandler.h>
#include<Common/TCostFunctionProvider.h>

template<typename Vc>
class TCouplingHandlerDensePrototype {
public:
	int xres,yres;
	int total;
	Vc *c;
	double *mu;
	TCouplingHandlerDensePrototype(int _xres, int _yres, Vc *_c, double *_mu);
	~TCouplingHandlerDensePrototype();
	inline int getRowLength(__attribute__((unused)) int x) {
		return yres;
	}
	inline int getColLength(__attribute__((unused)) int y) {
		return xres;
	}
	inline Vc getCRow(int x, int yIndex) {
		return c[yres*x+yIndex];
	}
	inline double getMuRow(int x, int yIndex) {
		return mu[yres*x+yIndex];
	}
	inline double getMuCol(int y, int xIndex) {
		return mu[yres*xIndex+y];
	}
	inline void setMuRow(int x, int yIndex, double value) {
		mu[yres*x+yIndex]=value;
	}
	void setMuRow(int x, double *valueList) {
		int i;
		for(i=0;i<getRowLength(x);i++) {
			mu[yres*x+i]=valueList[i];
		}
	}
	inline void setMuCol(int y, int xIndex, double value) {
		mu[yres*xIndex+y]=value;
	}
	inline int getColNr(__attribute__((unused)) int x, int yIndex) {
		return yIndex;
	}
	inline int getRowNr(__attribute__((unused)) int y, int xIndex) {
		return xIndex;
	}
	void clearMuRow(int x);
};


template<typename Vc>
class TCouplingHandlerSemiDensePrototype {
public:
	int xres,yres;
	int total;
	double *mu;
	Vc *c;
	TVarListHandler *xVars,*yVars;
	TCouplingHandlerSemiDensePrototype(int _xres, int _yres, Vc *_c, double *_mu, TVarListHandler *_xVars);
	~TCouplingHandlerSemiDensePrototype();
	inline int getRowLength(int x) {
		return xVars->lenList[x];
	}
	inline int getColLength(int y) {
		return yVars->lenList[y];
	}
	inline Vc getCRow(int x, int yIndex) {
		return c[yres*x+xVars->varList[x][yIndex]];
	}
	inline double getMuRow(int x, int yIndex) {
		return mu[yres*x+xVars->varList[x][yIndex]];
	}
	inline double getMuCol(int y, int xIndex) {
		return mu[yres*yVars->varList[y][xIndex]+y];
	}
	inline void setMuRow(int x, int yIndex, double value) {
		mu[yres*x+xVars->varList[x][yIndex]]=value;
	}
	void setMuRow(int x, double *valueList) {
		int i;
		for(i=0;i<getRowLength(x);i++) {
			mu[yres*x+xVars->varList[x][i] ]=valueList[i];
		}
	}
	inline void setMuCol(int y, int xIndex, double value) {
		mu[yres*yVars->varList[y][xIndex]+y]=value;
	}
	inline int getColNr(int x, int yIndex) {
		return xVars->varList[x][yIndex];
	}
	inline int getRowNr(int y, int xIndex) {
		return yVars->varList[y][xIndex];
	}
	
	void clearMuRow(int x);
	void updateXVars(TVarListHandler *_newXVars, bool keepXVars);

};




class TCouplingHandlerSparse {
public:
	int xres,yres;
	int total;
	double *mu;
	TCostFunctionProviderBase *cProvider;
	double *c;
	TVarListHandler *xVars;
	int *offsets;
	bool free_c;
	TCouplingHandlerSparse(int _xres, int _yres, TCostFunctionProviderBase *_cProvider, TVarListHandler *_xVars);

	~TCouplingHandlerSparse();
	
	void computeOffsets();
	
	inline int getRowLength(int x) {
		return xVars->lenList[x];
	}
	inline double getCRow(int x, int yIndex) {
		return c[offsets[x]+yIndex];
	}
	inline double getMuRow(int x, int yIndex) {
		return mu[offsets[x]+yIndex];
	}
	inline void setMuRow(int x, int yIndex, double value) {
		mu[offsets[x]+yIndex]=value;
	}
	void setMuRow(int x, double *valueList) {
		int i;
		for(i=0;i<getRowLength(x);i++) {
			mu[offsets[x]+i ]=valueList[i];
		}
	}
	inline int getColNr(int x, int yIndex) {
		return xVars->varList[x][yIndex];
	}
	
	void clearMuRow(int x);
	void updateXVars(TVarListHandler *_newXVars, bool keepXVars);
	
	double checkMarginalConstraints(double *muX, double *muY);

};


typedef TCouplingHandlerDensePrototype<double> TCouplingHandlerDense;
typedef TCouplingHandlerDensePrototype<int> TCouplingHandlerDenseInt;

typedef TCouplingHandlerSemiDensePrototype<double> TCouplingHandlerSemiDense;
typedef TCouplingHandlerSemiDensePrototype<int> TCouplingHandlerSemiDenseInt;


#endif
