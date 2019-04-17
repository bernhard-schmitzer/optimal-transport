#include<cstdlib>
#include <vector>

#include<Common.h>
#include<ShortCutSolver.h>
#include<LP_CPLEX.h>

int main() {
	int msg;
	
	// read raw data from file
	char filenamePosX[]="data/reflector_sphere/posX_9771.dat";
	char filenamePosY[]="data/reflector_sphere/posY_10800_Monge.dat";
	char filenameMuX[]="data/reflector_sphere/muX_9771.dat";
	char filenameMuY[]="data/reflector_sphere/muY_10800_Monge.dat";

	std::vector<double> posXdat=readFile<double>(filenamePosX);
	std::vector<double> posYdat=readFile<double>(filenamePosY);
	std::vector<double> muXdat=readFile<double>(filenameMuX);
	std::vector<double> muYdat=readFile<double>(filenameMuY);

	int dim=3;
	int xres=muXdat.size();
	int yres=muYdat.size();


	///////////////////////////////////////////////

	// setup problem data
	int posXdim[]={xres, dim};
	int posYdim[]={yres, dim};
	
	TDoubleMatrix posX;
	posX.data=posXdat.data();
	posX.dimensions=posXdim;
	posX.depth=2;

	TDoubleMatrix posY;
	posY.data=posYdat.data();
	posY.dimensions=posYdim;
	posY.depth=2;

	
	int depth=6;
	int nNeighbours=5;
	
	
	///////////////////////////////////////////////

	TMultiScaleSetupSingleSphere MultiScaleSetupX(&posX,muXdat.data(),depth);
	MultiScaleSetupX.HierarchyBuilderChildMode=THierarchyBuilder::CM_Tree;
	msg=MultiScaleSetupX.Setup();	
	if(msg!=0) {
		eprintf("%d\n",msg);
		return msg;
	}

	TMultiScaleSetupSingleSphere MultiScaleSetupY(&posY,muYdat.data(),depth);
	MultiScaleSetupY.HierarchyBuilderChildMode=THierarchyBuilder::CM_Tree;
	msg=MultiScaleSetupY.Setup();	
	if(msg!=0) {
		eprintf("%d\n",msg);
		return msg;
	}

	
	eprintf("hierarchical cardinalities:\n");
	for(int layer=0;layer<MultiScaleSetupX.nLayers;layer++) {
		eprintf("%d\t%d\n",layer,MultiScaleSetupX.HP->layers[layer]->nCells);
	}
	
	
	// first compute plain Euclidean radii, needed to determine nearest neighbours
	msg=MultiScaleSetupX.TMultiScaleSetupSingleBase::SetupRadii();	
	// compute neighbours for x, based on Euclidean distances
	MultiScaleSetupX.neighboursH=THierarchicalNN::getNeighboursH(
			MultiScaleSetupX.posH, MultiScaleSetupX.radii,
			MultiScaleSetupX.HP, nNeighbours);
	// free Euclidean radii data
	MultiScaleSetupX.HP->signal_free_double(MultiScaleSetupX.radii, 0, MultiScaleSetupX.nLayers-2);
	
	
	// project points back to sphere
	MultiScaleSetupX.SetupProjectPoints();
	MultiScaleSetupY.SetupProjectPoints();
	// compute sphere radii
	MultiScaleSetupX.SetupRadii();
	MultiScaleSetupY.SetupRadii();
	
	
	
	TCostFunctionProvider_Reflector costFunctionProvider(
			MultiScaleSetupX.resH, MultiScaleSetupY.resH,
			MultiScaleSetupX.posH, MultiScaleSetupY.posH,
			MultiScaleSetupX.nLayers, MultiScaleSetupX.dim
			);
	
			
	TShieldGeneratorTree_Reflector shieldGenerator(
			MultiScaleSetupX.dim,
			MultiScaleSetupY.HP, MultiScaleSetupY.posH, MultiScaleSetupY.radii,
			0 /* finest layer */, 0 /* coarsest layer */,
			MultiScaleSetupX.resH, MultiScaleSetupX.posH, MultiScaleSetupX.neighboursH
			);
	
	TShortCutCouplingHandlerInterface couplingHandlerInterface(
			MultiScaleSetupX.resH, MultiScaleSetupY.resH,
			MultiScaleSetupX.nLayers);
	
	
	// CPLEX
	TShortCutSubSolverInterfaceCPLEX subSolverInterface(
			MultiScaleSetupX.nLayers,
			MultiScaleSetupX.muH, MultiScaleSetupY.muH,
			&couplingHandlerInterface,true);
	
	
	///////////////////////////////////////////////

	TShortCutMultiScaleSolver MultiScaleSolver(
			&costFunctionProvider,
			&couplingHandlerInterface,
			&subSolverInterface,
			&shieldGenerator,
			MultiScaleSetupX.HP, MultiScaleSetupY.HP,
			1, // coarsest layer
			//TShortCutSolver::VCHECK_PRIMAL // need primal to test basis extraction
			TShortCutSolver::VCHECK_DUAL
			);
			
	MultiScaleSolver.autoDeletePointers=false;

	msg=MultiScaleSolver.solve();
	if(msg!=0) {
		eprintf("%d\n",msg);
		return msg;
	}
	
	
	return 0;

}
