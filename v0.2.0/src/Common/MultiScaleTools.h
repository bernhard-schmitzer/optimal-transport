#ifndef MultiScaleTools_H_
#define MultiScaleTools_H_

#include<vector>

#include<Common/PythonTypes.h>
#include<Common/ErrorCodes.h>
#include<Common/Tools.h>
#include<Common/Verbose.h>

#include<Common/GridTools.h>
#include<Common/THierarchyBuilder.h>
#include<Common/THierarchicalPartition.h>
#include<Common/THierarchicalCostFunctionProvider.h>


//class TMultiScaleSetupBase {
//public:
//	TDoubleMatrix *posX, *posY; // point clouds for marginal positions
//	double *muX, *muY; // pointers to marginal masses
//	int xres, yres; // integers for total cardinality of marginal points
//	int dim; // dimensionality of marginal positions
//	
//	// hierarchy
//	int depth; // parameter for controlling nr of layers in hierarchical partition
//	int nLayers; // number of layers (usually depth+2)
//	THierarchyBuilder *HBX,*HBY; // hierarchy builder classes
//	THierarchicalPartition *HPX,*HPY; // hierarchical partition classes
//	double **posXH, **posYH; // hierarchical positions
//	double **muXH,**muYH; // hierarchical masses
//	int *xresH, *yresH; // hierarchical marginal cardinality
//	int HierarchyBuilderChildMode; // mode for handling children nodes in hierarchy builder
//		// set to THierarchyBuilder::CM_Grid in constructor.
//		
//	// other hierarchical stuff (setup on demand)
//	double **alphaH, **betaH; // dual potentials
//	double **xRadii, **yRadii; // radii for hiearchical partition cells
//	// list of neighbour points for X
//	TVarListHandler **xNeighboursH;
//		
//	//////////////////////////////////////////////////////////////////////////////
//	
//	TMultiScaleSetupBase(TDoubleMatrix *_posX, TDoubleMatrix *_posY, double *_muX, double *_muY, int _depth);
//	virtual ~TMultiScaleSetupBase();

//	int BasicSetup();
//	int BasicMeasureChecks();
//	int SetupHierarchicalPartition(double *mu, double *pos, int res, int dim, int depth,
//			THierarchyBuilder **_HB, THierarchicalPartition **_HP, double ***_posH, double ***_muH, int **_resH);
//	virtual int SetupHierarchicalPartitions();
//	
//	virtual int Setup();
//	virtual int SetupDuals();
//	virtual int SetupRadii();
//	
//};


//class TMultiScaleSetupGrid : public TMultiScaleSetupBase {
//public:
//	// multidimensional arrays of marginal measures
//	TDoubleMatrix *muXGrid, *muYGrid;
//	// grid dimensions of each hierarchy level. required for shield generators
//	// these are contiguous flattened 2d arrays with dimensions nLayers*dim
//	int *xDimH,*yDimH;
//	
//	TMultiScaleSetupGrid(TDoubleMatrix *_muXGrid, TDoubleMatrix *_muYGrid, int _depth);
//	int SetupHierarchicalPartitions();
//	int SetupGridNeighboursX();
//	// destructor
//	virtual ~TMultiScaleSetupGrid();
//	
//};

// Barycenter Version

class TMultiScaleSetupBarycenterBase {
public:
	int nMarginals; // number of marginals
	TDoubleMatrix **pos, *posZ; // point clouds for marginal positions and for center
	double **mu, *muZ; // pointers to marginal masses and central reference measure
	int *res; // cardinality of marginals
	int zres; // cardinality of center
	int dim; // dimensionality of positions
	
	// hierarchy
	int depth; // parameter for controlling nr of layers in hierarchical partition
	int nLayers; // number of layers (usually depth+2)
	THierarchyBuilder **HB,*HBZ; // hierarchy builder classes
	THierarchicalPartition **HP,*HPZ; // hierarchical partition classes
	double ***posH, **posZH; // hierarchical positions
	double ***muH,**muZH; // hierarchical masses
	int **resH, *zresH; // hierarchical marginal cardinality
	int HierarchyBuilderChildMode; // mode for handling children nodes in hierarchy builder
		// set to THierarchyBuilder::CM_Grid in constructor.
		
	// other hierarchical stuff (setup on demand)
	double ***alphaH, ***betaH; // dual potentials
	double ***radii, **zRadii; // radii for hiearchical partition cells
		
	//////////////////////////////////////////////////////////////////////////////
	
	TMultiScaleSetupBarycenterBase(
			int _nMarginals,
			TDoubleMatrix **_pos, TDoubleMatrix *_posZ,
			double **_mu, double *_muZ,
			int _depth);
	virtual ~TMultiScaleSetupBarycenterBase();

	int BasicSetup();
	int BasicMeasureChecks();
	int SetupHierarchicalPartition(double *aMu, double *aPos, int aRes, int dim, int depth,
			THierarchyBuilder **_aHB, THierarchicalPartition **_aHP, double ***_aPosH, double ***_aMuH, int **_aResH);
	int SetupHierarchicalPartitions();
	
	virtual int Setup();
	virtual int SetupDuals();
	virtual int SetupRadii();
	
};

// Algorithm for hierarchical nearest neighbour search

class THierarchicalNN {
public:
	struct TCandidate {
		int layer; // on which layer
		int z; // row/column index
		double dist; // lower bound on distance
		bool operator<(const THierarchicalNN::TCandidate& rhs) const { return dist < rhs.dist; };
	};
	
	typedef THierarchicalSearchList<THierarchicalNN::TCandidate> TCandidateList;


	static std::vector<int> find(double *posX, double **posYH, double **radiiY,
			THierarchicalPartition *HPY, int layerBottom, int nElements);
		// find nElement nearest neighbours of posX in the points of HPY
		// (with coordinates posYH and cell radii radiiY) at layer layerBottom
	
	static TVarListHandler* getNeighbours(double **posXH, double **radiiX,
			THierarchicalPartition *HPX, int layerBottom, int nElements);
		// find nElement nearest neighbours for each point in layerBottom of HPX
		// (with coordinates posXH, cell radii radiiX) by using the above find(...)
		// function on each point

	static TVarListHandler** getNeighboursH(double **posXH, double **radiiX,
			THierarchicalPartition *HPX, int nElements);
		// find nElement nearest neighbours for each point in layerBottom of HPX
		// (with coordinates posXH, cell radii radiiX) by using the above find(...)
		// function on each point

};


class THierarchicalDualMaximizer {
public:
	static constexpr int MODE_ALPHA=0;
	static constexpr int MODE_BETA=1;
	
	struct TCandidate {
		int layer; // on which layer
		int z; // row/column index
		double v; // uper bound on dual variable
		bool operator<(const THierarchicalDualMaximizer::TCandidate& rhs) const { return v < rhs.v; };
	};

	typedef THierarchicalSearchList<THierarchicalDualMaximizer::TCandidate> TCandidateList;

	static void getMaxDual(THierarchicalPartition *partitionX, THierarchicalPartition *partitionY,
			double **alpha, double **beta, int layerFine,
			THierarchicalCostFunctionProvider *costProvider,
			int mode);


};


// Single marginal setup

class TMultiScaleSetupSingleBase {
public:
	TDoubleMatrix *pos; // point clouds for marginal positions
	double *mu; // pointers to marginal masses
	int res; // integers for total cardinality of marginal points
	int dim; // dimensionality of marginal positions
	
	// hierarchy
	int depth; // parameter for controlling nr of layers in hierarchical partition
	int nLayers; // number of layers (usually depth+2)
	THierarchyBuilder *HB; // hierarchy builder classes
	THierarchicalPartition *HP; // hierarchical partition classes
	double **posH; // hierarchical positions
	double **muH; // hierarchical masses
	int *resH; // hierarchical marginal cardinality
	int HierarchyBuilderChildMode; // mode for handling children nodes in hierarchy builder
		// set to THierarchyBuilder::CM_Grid in constructor.
		
	// other hierarchical stuff (setup on demand)
	double **alphaH; // dual potentials
	double **radii; // radii for hiearchical partition cells
	// list of neighbour points for X
	TVarListHandler **neighboursH;
		
	//////////////////////////////////////////////////////////////////////////////
	
	TMultiScaleSetupSingleBase(TDoubleMatrix *_pos, double *_mu, int _depth);
	virtual ~TMultiScaleSetupSingleBase();

	int BasicSetup();
	int BasicMeasureChecks();
	virtual int SetupHierarchicalPartition();
	
	virtual int Setup();
	virtual int SetupDuals();
	virtual int SetupRadii();
	
};

// Cartesian grid version

class TMultiScaleSetupSingleGrid : public TMultiScaleSetupSingleBase {
public:
	// multidimensional array of marginal measure
	TDoubleMatrix *muGrid;
	// grid dimensions of each hierarchy level. required for shield generators
	// this is a contiguous flattened 2d array with dimensions nLayers*dim
	int *dimH;
	
	TMultiScaleSetupSingleGrid(TDoubleMatrix *_muGrid, int _depth);
	int SetupHierarchicalPartition();
	int SetupGridNeighbours();
	// destructor
	virtual ~TMultiScaleSetupSingleGrid();
	
};


class TMultiScaleSetupBarycenterContainer {
// class that provides arrays for relevant fields of several TMultiScaleSetupSingleBase instances
// to use these in a SinkhornBarycenter solver
public:
	int nMarginals;
	THierarchicalPartition **HP,*HPZ; // hierarchical partition classes
	double ***muH,**muZH; // hierarchical masses
	double ***alphaH, ***betaH; // dual potentials
	int **resH, *resZH; // cardinalities of layers
	double *weights;
	THierarchicalCostFunctionProvider **costProvider;

	TMultiScaleSetupBarycenterContainer();
	TMultiScaleSetupBarycenterContainer(const int _nMarginals);
	~TMultiScaleSetupBarycenterContainer();
	void setupEmpty(const int _nMarginals);
	void cleanup();
	void setMarginal(const int n, TMultiScaleSetupSingleBase &multiScaleSetup, const double weight);
	void setCenterMarginal(TMultiScaleSetupSingleBase &multiScaleSetup);
	void setCostFunctionProvider(const int n, THierarchicalCostFunctionProvider &costFunctionProvider);
};

#endif

