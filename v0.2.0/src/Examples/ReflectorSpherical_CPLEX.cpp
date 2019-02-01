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

	
	int depth=5;
	int nNeighbours=5;
	
	
	///////////////////////////////////////////////

	TMultiScaleSetupSphere MultiScaleSetup(&posX,&posY,muXdat.data(),muYdat.data(),depth);
	MultiScaleSetup.HierarchyBuilderChildMode=THierarchyBuilder::CM_Tree;
	msg=MultiScaleSetup.Setup();	
	if(msg!=0) {
		eprintf("%d\n",msg);
		return msg;
	}
	
	eprintf("hierarchical cardinalities:\n");
	for(int layer=0;layer<MultiScaleSetup.nLayers;layer++) {
		eprintf("%d\t%d\n",layer,MultiScaleSetup.HPX->layers[layer]->nCells);
	}
	
	
	// first compute plain Euclidean radii, needed to determine nearest neighbours
	msg=MultiScaleSetup.TMultiScaleSetupBase::SetupRadii();	
	// compute neighbours for x, based on Euclidean distances
	MultiScaleSetup.xNeighboursH=THierarchicalNN::getNeighboursH(
			MultiScaleSetup.posXH, MultiScaleSetup.xRadii,
			MultiScaleSetup.HPX, nNeighbours);
	// free Euclidean radii data
	MultiScaleSetup.HPX->signal_free_double(MultiScaleSetup.xRadii, 0, MultiScaleSetup.nLayers-2);
	MultiScaleSetup.HPY->signal_free_double(MultiScaleSetup.yRadii, 0, MultiScaleSetup.nLayers-2);
	
	
	// project points back to sphere
	MultiScaleSetup.SetupProjectPoints();
	// compute sphere radii
	MultiScaleSetup.SetupRadii();
	
	
	
	TCostFunctionProvider_Reflector costFunctionProvider(
			MultiScaleSetup.xresH, MultiScaleSetup.yresH,
			MultiScaleSetup.posXH, MultiScaleSetup.posYH,
			MultiScaleSetup.nLayers, MultiScaleSetup.dim
			);
	
			
	TShieldGeneratorTree_Reflector shieldGenerator(
			MultiScaleSetup.dim,
			MultiScaleSetup.HPY, MultiScaleSetup.posYH, MultiScaleSetup.yRadii,
			0 /* finest layer */, 0 /* coarsest layer */,
			MultiScaleSetup.xresH, MultiScaleSetup.posXH, MultiScaleSetup.xNeighboursH
			);
	
	TShortCutCouplingHandlerInterface couplingHandlerInterface(
			MultiScaleSetup.xresH, MultiScaleSetup.yresH,
			MultiScaleSetup.nLayers);
	
	
	// CPLEX
	TShortCutSubSolverInterfaceCPLEX subSolverInterface(
			MultiScaleSetup.nLayers,
			MultiScaleSetup.muXH, MultiScaleSetup.muYH,
			&couplingHandlerInterface,true);
	
	
	///////////////////////////////////////////////

	TShortCutMultiScaleSolver MultiScaleSolver(
			&costFunctionProvider,
			&couplingHandlerInterface,
			&subSolverInterface,
			&shieldGenerator,
			MultiScaleSetup.HPX, MultiScaleSetup.HPY,
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
