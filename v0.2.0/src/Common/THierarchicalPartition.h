#ifndef THierarchicalPartition_H_
#define THierarchicalPartition_H_

#include<cstdlib>
#include"Common/TVarListHandler.h"

class TPartitionLayer {
public:
	int nCells; // number of cells in that layer
	// only the parent-children-leaves structure is fundamental. the rest can be encoded as "signals" living on the partition
	// implement these signals separately
	int *parent;
	int **children, **leaves;
	int *nChildren, *nLeaves;
	TPartitionLayer();
	~TPartitionLayer();
	void initializeEmpty(int _nCells);
};

class THierarchicalPartition {
public:
	static const int MODE_MIN=0;
	static const int MODE_MAX=1;
	TPartitionLayer **layers;
	int nLayers; // number of layers
	int dim; // dimension of ambient space
	THierarchicalPartition(int _nLayers, int _dim);
	~THierarchicalPartition();
	
	void computeHierarchicalMasses(double *mu, double **muLayers);
	double** signal_allocate_double(int lTop, int lBottom);
	void signal_free_double(double **signal, int lTop, int lBottom);
	void signal_propagate_double(double **signal ,int lTop, int lBottom, int mode);
	
	void signal_refine_double(double *signal, double *signalFine, int lTop);
};


template<typename T>
class THierarchicalProductSignal {
public:
	static const int MODE_MIN=0;
	static const int MODE_MAX=1;
	
	// temporariy "global" variables for check_dualConstraints
	T **c, **alpha, **beta;
	T slack;
	TVarListHandler* varList;
	// for advanced constraint checking (possibly need to clean this up later)
	T **slackOffsetX, **slackOffsetY;
	
	THierarchicalPartition *partitionX, *partitionY;
	THierarchicalProductSignal(THierarchicalPartition *_partitionX, THierarchicalPartition *_partitionY);
	void signal_propagate(T **signal, int lTop, int lBottom, int mode);
	TVarListHandler* check_dualConstraints(T **_c, T **_alpha, T **_beta, int lTop, int lBottom, T _slack);
	void check_dualConstraints_iterateTile(int l, int x, int y, int lBottom);
	//
	TVarListHandler* check_dualConstraints_adaptive(T **_c, T **_alpha, T **_beta, int lTop, int lBottom,
		T **_slackOffsetX, T **_slackOffsetY);
	void check_dualConstraints_adaptive_iterateTile(int l, int x, int y, int lBottom);
};


TVarListHandler* refineVarList(THierarchicalPartition *partitionX, THierarchicalPartition *partitionY,
		TVarListHandler *varListCoarse, int layerIdCoarse, bool doSort=false);

#endif
