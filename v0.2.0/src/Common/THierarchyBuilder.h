#ifndef THierarchyBuilder_H_
#define THierarchyBuilder_H_

#include<cstdlib>
#include<cmath>
#include<vector>
#include<Common/THierarchicalPartition.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class THierarchyBuilderNode {
public:
	std::vector<int> children;
	std::vector<int> leaves;
	std::vector<int> posCode;
	int parent;
};

class THierarchyBuilderLayer {
public:
	std::vector<THierarchyBuilderNode> nodes;
};	


class THierarchyBuilder {
public:
	static constexpr int CM_Tree=0;
	static constexpr int CM_Grid=1;
	static constexpr double boxTolerance=1.E-10;
	
	double *points;
	int nPoints; // number of points at lowest level
	int dim;
	std::vector<double> boxLo, boxHi; // lower and higher bounds for box around point cloud;
	std::vector<THierarchyBuilderLayer> layers;
	int childMode;
	
	THierarchyBuilder(double *_points, int _nPoints, int _dim,
			int _childMode, int partitionDepth);
	
	void setBox();
	void setRoot();
	void refine();
	std::vector<std::vector<int> > getChildrenPosCodes(int layerId, int nodeId);
	std::vector<std::vector<int> > getChildrenLeaves(int layerId, int nodeId,
		const std::vector<std::vector<int> >& childrenPosCodes);
	void getRelPosCodeFromIndex(int index, int dim, int *posCode);
	void getOffsetPosCode(int *relPosCode, int *parentPosCode, int dim);
	bool isInBox(double *coord, const int * const posCode, int dim, int layerId);
	void addAtomicLayer();
	
	static double max(double *x, int n, int step, int offset);
	static double min(double *x, int n, int step, int offset);
	
	THierarchicalPartition* convert();
	double** allocateDoubleSignal(int sigdim, int nLayers=0);
	void freeSignal(double **signal, int nLayers);
	void getSignalPos(double **signal);
	int* getDimH(int *finestDims);
	int* getResH();
	
	double** getSignalRadii();
	double** getSignalRadiiAxis(int *axis, int naxis);
};


#endif
