#include<cstdlib>
#include <vector>

#include<Common.h>
#include<Sinkhorn.h>

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
	int layerCoarsest=0;
	int layerFinest=depth;

	
	///////////////////////////////////////////////

	TMultiScaleSetupSingleSphere MultiScaleSetupX(&posX,muXdat.data(),depth);
	MultiScaleSetupX.HierarchyBuilderChildMode=THierarchyBuilder::CM_Tree;
	msg=MultiScaleSetupX.Setup();	
	if(msg!=0) { eprintf("error: %d\n",msg); return msg; }	
	// project points back to sphere
	MultiScaleSetupX.SetupProjectPoints();
	// compute sphere radii
	msg=MultiScaleSetupX.SetupRadii();
	if(msg!=0) { eprintf("error: %d\n",msg); return msg; }
	// allocate hierarchical dual variables
	msg=MultiScaleSetupX.SetupDuals();
	if(msg!=0) { eprintf("error: %d\n",msg); return msg; }


	TMultiScaleSetupSingleSphere MultiScaleSetupY(&posY,muYdat.data(),depth);
	MultiScaleSetupY.HierarchyBuilderChildMode=THierarchyBuilder::CM_Tree;
	msg=MultiScaleSetupY.Setup();	
	if(msg!=0) { eprintf("error: %d\n",msg); return msg; }	
	// project points back to sphere
	MultiScaleSetupY.SetupProjectPoints();
	// compute sphere radii
	msg=MultiScaleSetupY.SetupRadii();
	if(msg!=0) { eprintf("error: %d\n",msg); return msg; }
	// allocate hierarchical dual variables
	msg=MultiScaleSetupY.SetupDuals();
	if(msg!=0) { eprintf("error: %d\n",msg); return msg; }


	
	eprintf("hierarchical cardinalities:\n");
	for(int layer=0;layer<MultiScaleSetupX.nLayers;layer++) {
		eprintf("%d\t%d\n",layer,MultiScaleSetupX.HP->layers[layer]->nCells);
	}
	
	
	/////////////////////////////////////////////////////////////////////////////////////////
	// setup a cost function provider
	THierarchicalCostFunctionProvider_Reflector costProvider(
		MultiScaleSetupX.posH, MultiScaleSetupY.posH,
		MultiScaleSetupX.radii, MultiScaleSetupY.radii,
		dim, 0,
		true,
		MultiScaleSetupX.alphaH, MultiScaleSetupY.alphaH
		);
	
	
	/////////////////////////////////////////////////////////////////////////////////////////
	// epsScaling
	double epsStart=1E1;
	double epsTarget=2E-5;
	int epsSteps=25;
	double epsBoxScale=.7;
	
	TEpsScalingHandler epsScalingHandler(epsStart,epsTarget,epsSteps); // basic eps scaling
	epsScalingHandler.getEpsScalesFromBox(epsBoxScale,2.,MultiScaleSetupX.nLayers); // eps scales for each layer
	msg=epsScalingHandler.getEpsScalingSplit(layerCoarsest,1); // create sub eps lists
	if(msg!=0) { eprintf("error: %d\n",msg); return msg; }

	
	/////////////////////////////////////////////////////////////////////////////////////////
	// other parameters
	TSinkhornSolverBase::TSinkhornSolverParameters cfg={
			1E-5, // maxError
			100000, // maxIterations
			100, // innerIterations
			20, // maxAbsorptionLoops
			1E3, // absorption_scalingBound
			1E3, // absorption_scalingLowerBound
			1E-6, // truncation_thresh
			//10., // turncation_thresh
			true // refineKernel
			};

	/////////////////////////////////////////////////////////////////////////////////////////
	// create solver object
	TSinkhornSolverStandard SinkhornSolver(MultiScaleSetupX.nLayers, epsScalingHandler.nEpsLists, epsScalingHandler.epsLists,
			layerCoarsest, layerFinest,
			cfg,
			MultiScaleSetupX.HP, MultiScaleSetupY.HP,
			MultiScaleSetupX.muH, MultiScaleSetupY.muH,
			MultiScaleSetupX.alphaH, MultiScaleSetupY.alphaH,
			&costProvider
			);
	

	SinkhornSolver.initialize();
	msg=SinkhornSolver.solve();	
	printf("return code: %d\n",msg);

	
	// recompute kernel one last time
	SinkhornSolver.generateKernel();
	// return primal objective value
	double primalScore=SinkhornSolver.scorePrimalUnreg();
	printf("primal score: %e\n",primalScore);

	
	return 0;

}
