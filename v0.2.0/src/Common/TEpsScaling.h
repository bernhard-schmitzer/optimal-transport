#ifndef TEPSSCALING_H_
#define TEPSSCALING_H_

#include<cstdlib>
#include<vector>
#include<cmath>
#include<cstring>

#include<Common/ErrorCodes.h>
#include<Common/Verbose.h>


class TEpsScalingHandler {
public:
	// basic eps scaling: largest, smallest eps, number of steps for decrease
	double epsStart, epsTarget;
	int epsSteps;
	double *epsList; // list of all eps values
	
	// hierarchical eps scaling over multiple layers
	double *epsScales; // lower bound for eps at each layer
	int nLayers; // number of layers
	double **epsLists; // one eps list per layer
	int *nEpsLists; // number of eps in each layer

	
	
	TEpsScalingHandler(double _epsStart, double _epsTarget, int _epsSteps);
	~TEpsScalingHandler();
	void getEpsScalesFromBox(double boxScale, double layerExponent, int _nLayers);
	int getEpsScalingSplit(int levelCoarsest, int nOverlap);

	void getEpsScalingAllFinest(int levelFinest);
};

#endif
