#include<cstdlib>
#include <vector>

#include<Common.h>
#include<ShortCutSolver.h>

#ifdef USE_CPLEX
#include<LP_CPLEX.h>
#endif
#ifdef USE_LEMON
#include<LP_Lemon.h>
#endif
#ifdef USE_LPSOLVE
#include<LP_lp_solve.h>
#endif


int main(int argc, char* argv[]);
int W2Grid();
int WpGrid();
int W2Grid_64x64();
int WpSphere();



int main(int argc, char* argv[]) {
	char option;
	if(argc<=1) {
		option='0';
	} else {
		option=argv[1][0];
	}
	switch(option) {
		case '1':
			WpGrid();
			break;
		case '2':
			W2Grid_64x64();
			break;
		case '3':
			WpSphere();
			break;
		default:
			W2Grid();
			break;
	}			
}


int W2Grid() {
	int msg;
	
	// both marginals life on a regular 8x8 grid
	// provide measure data as consecutive 64 double entries
	double muXdat[]={6.101275e-03, 7.204110e-03, 8.535062e-03, 9.910280e-03, 1.104128e-02, 1.166015e-02, 1.167389e-02, 1.123175e-02, 6.704433e-03, 8.128397e-03, 9.854252e-03, 1.163492e-02, 1.309322e-02, 1.388739e-02, 1.390795e-02, 1.335771e-02, 7.292185e-03, 9.029484e-03, 1.114963e-02, 1.334184e-02, 1.513976e-02, 1.613150e-02, 1.620053e-02, 1.561280e-02, 7.908921e-03, 9.950044e-03, 1.245024e-02, 1.503625e-02, 1.715816e-02, 1.834562e-02, 1.848901e-02, 1.792387e-02, 8.567936e-03, 1.091527e-02, 1.379087e-02, 1.675992e-02, 1.919501e-02, 2.057699e-02, 2.080940e-02, 2.028992e-02, 9.247266e-03, 1.191367e-02, 1.517757e-02, 1.854434e-02, 2.131180e-02, 2.291030e-02, 2.324255e-02, 2.273922e-02, 9.914861e-03, 1.291900e-02, 1.659968e-02, 2.040437e-02, 2.355167e-02, 2.540296e-02, 2.582486e-02, 2.523923e-02, 1.055807e-02, 1.392742e-02, 1.806860e-02, 2.236741e-02, 2.594533e-02, 2.805923e-02, 2.849642e-02, 2.764292e-02};

	double muYdat[]={3.908070e-03, 4.687770e-03, 5.569132e-03, 6.429051e-03, 7.125903e-03, 7.561235e-03, 7.732784e-03, 7.742170e-03, 5.244827e-03, 5.949688e-03, 6.792079e-03, 7.684573e-03, 8.490201e-03, 9.083672e-03, 9.421825e-03, 9.576329e-03, 7.285414e-03, 7.776036e-03, 8.404373e-03, 9.162126e-03, 9.953290e-03, 1.064454e-02, 1.115002e-02, 1.149372e-02, 1.034929e-02, 1.044027e-02, 1.058721e-02, 1.093244e-02, 1.148517e-02, 1.214090e-02, 1.277317e-02, 1.333126e-02, 1.502886e-02, 1.456200e-02, 1.391201e-02, 1.346945e-02, 1.344750e-02, 1.383385e-02, 1.447324e-02, 1.520699e-02, 2.202719e-02, 2.099936e-02, 1.930211e-02, 1.767088e-02, 1.664565e-02, 1.640100e-02, 1.678885e-02, 1.752175e-02, 3.157802e-02, 3.026055e-02, 2.748538e-02, 2.437437e-02, 2.190360e-02, 2.055179e-02, 2.025589e-02, 2.062372e-02, 4.272636e-02, 4.166590e-02, 3.810900e-02, 3.349480e-02, 2.927364e-02, 2.633617e-02, 2.482542e-02, 2.433618e-02};
	
	// dimensions of grid
	int muXdim[]={8, 8};
	int muYdim[]={8, 8};
	
	
	// put raw marginal data into small container structs
	TDoubleMatrix muX,muY;

	muX.data=muXdat;
	muX.dimensions=muXdim;
	muX.depth=2;
	
	muY.data=muYdat;
	muY.dimensions=muYdim;
	muY.depth=2;
	
	// fundamental parameters
	int depth=3; // hierarchy depth

	#ifdef USE_LEMON
	// Lemon parameters	
	double measureScale=1E-9; // scale to which marginal measures are truncated
	double cScale=1E-7; // scale to which cost function is truncated
	bool dualOffset=true; // activate small trick with dual variables
		// that slightly accelerates the lemon network simplex
	
	///////////////////////////////////////////////
	// total nr of points in each marginal
	int xres=GridToolsGetTotalPoints(muX.depth, muX.dimensions);
	int yres=GridToolsGetTotalPoints(muY.depth, muY.dimensions);

	// truncate marginals
	// both marginals are rounded to multiples of measureScale
	msg=MeasureToolsTruncateMeasures(muX.data, muY.data, xres, yres, measureScale);
	if(msg!=0) {
		return msg;
	}
	#endif


	///////////////////////////////////////////////
	// problem setup
	
	
	// generate multi-scale problem representation:
	// the MultiScaleSetup object will subsequently hold all
	// required data for the solver and other algorithms
	
	// the class TMultiScaleSetupGrid assumes that muX and muY describe measures
	// that live on regular Cartesian grids with edge lenghts 1
	
	// the more general base class TMultiScaleSetupBasefor uses marginals with support on arbitrary point clouds
	TMultiScaleSetupSingleGrid MultiScaleSetupX(&muX,depth);
	msg=MultiScaleSetupX.Setup();
	if(msg!=0) {
		printf("%d\n",msg);
		return msg;
	}


	TMultiScaleSetupSingleGrid MultiScaleSetupY(&muY,depth);
	msg=MultiScaleSetupY.Setup();
	if(msg!=0) {
		printf("%d\n",msg);
		return msg;
	}


	
	// setup various aux components for actual solver algorithm
	TCostFunctionProvider_Dynamic costFunctionProvider(
			MultiScaleSetupX.resH, MultiScaleSetupY.resH,
			MultiScaleSetupX.posH, MultiScaleSetupY.posH,
			MultiScaleSetupX.nLayers, MultiScaleSetupX.dim);
	
	TShieldGeneratorGrid_SqrEuclidean shieldGenerator(
			MultiScaleSetupX.dim,
			MultiScaleSetupX.dimH, MultiScaleSetupY.dimH,
			MultiScaleSetupX.nLayers);
			
	TShortCutCouplingHandlerInterface couplingHandlerInterface(
			MultiScaleSetupX.resH, MultiScaleSetupY.resH,
			MultiScaleSetupX.nLayers);
	
	
	#ifdef USE_CPLEX
	// LP subsolver: CPLEX 
	TShortCutSubSolverInterfaceCPLEX subSolverInterface(
			MultiScaleSetupX.nLayers,
			MultiScaleSetupX.muH, MultiScaleSetupY.muH,
			&couplingHandlerInterface,true);
	
	static constexpr int VIOLATION_CHECKMODE=TShortCutSolver::VCHECK_DUAL;		
	
	#endif
	#ifdef USE_LEMON
	// LP subsolver: Lemon
	TShortCutSubSolverInterfaceLemon subSolverInterface(
			MultiScaleSetupX.nLayers,
			MultiScaleSetupX.muH, MultiScaleSetupY.muH,
			&couplingHandlerInterface,
			measureScale, cScale,
			dualOffset);

	static constexpr int VIOLATION_CHECKMODE=TShortCutSolver::VCHECK_PRIMAL;

	#endif
	#ifdef USE_LPSOLVE
	// LP subsolver: lp_solve 
	TShortCutSubSolverInterfaceLpSolve subSolverInterface(
			MultiScaleSetupX.nLayers,
			MultiScaleSetupX.muH, MultiScaleSetupY.muH,
			&couplingHandlerInterface,true);
	
	static constexpr int VIOLATION_CHECKMODE=TShortCutSolver::VCHECK_DUAL;		
	#endif


	///////////////////////////////////////////////
	// main solver object
	TShortCutMultiScaleSolver MultiScaleSolver(
			&costFunctionProvider,
			&couplingHandlerInterface,
			&subSolverInterface,
			&shieldGenerator,
			MultiScaleSetupX.HP, MultiScaleSetupY.HP,
			1, // coarsest layer
			VIOLATION_CHECKMODE
			);
			
	MultiScaleSolver.autoDeletePointers=false;

	msg=MultiScaleSolver.solve();
	if(msg!=0) {
		printf("%d\n",msg);
		return msg;
	}
	
	// just demonstrate a few things that can be done after solving
	printf("\n");
	printf("Wasserstein-2 distance: %f\n",std::pow(MultiScaleSolver.objective,1./2.));
	
	printf("extract optimal coupling\n");
	// extract optimal coupling in sparse pos/coo format (see e.g. scipy documentation)
	TSparsePosContainer coupling=MultiScaleSolver.couplingHandlerInterface->getCouplingPos();
	for(int i=0;i<coupling.nParticles;i++) {
		printf("\t%d\t%d\t%f\n",coupling.posStart[i],coupling.posEnd[i],coupling.mass[i]);
	}
	
	printf("optimal dual variables\n");
	printf("alpha\n");
	for(int i=0;i<MultiScaleSetupX.res;i++) {
		printf("\t%f\n",MultiScaleSolver.alpha[i]);
	}
	printf("beta\n");
	for(int i=0;i<MultiScaleSetupY.res;i++) {
		printf("\t%f\n",MultiScaleSolver.beta[i]);
	}
	

	// verify dual constraints
	printf("Verify all dual constraints\n");
	
	double cost;
	bool violated=false;
	for(int x=0;x<MultiScaleSetupX.res;x++) {
		for(int y=0;y<MultiScaleSetupY.res;y++) {
			cost=costFunctionProvider.getCValue(x,y);
			#ifdef USE_LEMON
			cost=((int) (cost/cScale))*cScale;
			#endif
			if(MultiScaleSolver.alpha[x]+MultiScaleSolver.beta[y]>=cost+1E-13) {
				printf("%d\t%d\t%f\n",x,y,MultiScaleSolver.alpha[x]+MultiScaleSolver.beta[y]-cost);
				violated=true;
			}
		}
	}
	
	if(violated) {
		printf("not all dual constraints are satisfied! something may have gone wrong!\n");
		printf("(or it's just a truncation artifact of the lemon solver)\n");
	} else {
		printf("all dual constraints are satisfied\n");
	}

	
	return 0;

}




int WpGrid() {
	int msg;
	
	// set up marginals as in W2Grid
	double muXdat[]={6.101275e-03, 7.204110e-03, 8.535062e-03, 9.910280e-03, 1.104128e-02, 1.166015e-02, 1.167389e-02, 1.123175e-02, 6.704433e-03, 8.128397e-03, 9.854252e-03, 1.163492e-02, 1.309322e-02, 1.388739e-02, 1.390795e-02, 1.335771e-02, 7.292185e-03, 9.029484e-03, 1.114963e-02, 1.334184e-02, 1.513976e-02, 1.613150e-02, 1.620053e-02, 1.561280e-02, 7.908921e-03, 9.950044e-03, 1.245024e-02, 1.503625e-02, 1.715816e-02, 1.834562e-02, 1.848901e-02, 1.792387e-02, 8.567936e-03, 1.091527e-02, 1.379087e-02, 1.675992e-02, 1.919501e-02, 2.057699e-02, 2.080940e-02, 2.028992e-02, 9.247266e-03, 1.191367e-02, 1.517757e-02, 1.854434e-02, 2.131180e-02, 2.291030e-02, 2.324255e-02, 2.273922e-02, 9.914861e-03, 1.291900e-02, 1.659968e-02, 2.040437e-02, 2.355167e-02, 2.540296e-02, 2.582486e-02, 2.523923e-02, 1.055807e-02, 1.392742e-02, 1.806860e-02, 2.236741e-02, 2.594533e-02, 2.805923e-02, 2.849642e-02, 2.764292e-02};

	double muYdat[]={3.908070e-03, 4.687770e-03, 5.569132e-03, 6.429051e-03, 7.125903e-03, 7.561235e-03, 7.732784e-03, 7.742170e-03, 5.244827e-03, 5.949688e-03, 6.792079e-03, 7.684573e-03, 8.490201e-03, 9.083672e-03, 9.421825e-03, 9.576329e-03, 7.285414e-03, 7.776036e-03, 8.404373e-03, 9.162126e-03, 9.953290e-03, 1.064454e-02, 1.115002e-02, 1.149372e-02, 1.034929e-02, 1.044027e-02, 1.058721e-02, 1.093244e-02, 1.148517e-02, 1.214090e-02, 1.277317e-02, 1.333126e-02, 1.502886e-02, 1.456200e-02, 1.391201e-02, 1.346945e-02, 1.344750e-02, 1.383385e-02, 1.447324e-02, 1.520699e-02, 2.202719e-02, 2.099936e-02, 1.930211e-02, 1.767088e-02, 1.664565e-02, 1.640100e-02, 1.678885e-02, 1.752175e-02, 3.157802e-02, 3.026055e-02, 2.748538e-02, 2.437437e-02, 2.190360e-02, 2.055179e-02, 2.025589e-02, 2.062372e-02, 4.272636e-02, 4.166590e-02, 3.810900e-02, 3.349480e-02, 2.927364e-02, 2.633617e-02, 2.482542e-02, 2.433618e-02};

	int muXdim[]={8, 8};
	int muYdim[]={8, 8};
	
	TDoubleMatrix muX,muY;

	muX.data=muXdat;
	muX.dimensions=muXdim;
	muX.depth=2;
	
	muY.data=muYdat;
	muY.dimensions=muYdim;
	muY.depth=2;
	


	// fundamental parameters
	int depth=3; // hierarchy depth
	double p=1.5; // exponent for Wp distance that we want to compute
		// must be strictly >1., becomes slower as p approaches 1
	
	#ifdef USE_LEMON
	// Lemon parameters	
	double measureScale=1E-9; // scale to which marginal measures are truncated
	double cScale=1E-7; // scale to which cost function is truncated
	bool dualOffset=true; // activate small trick with dual variables
		// that slightly accelerates the lemon network simplex
	
	///////////////////////////////////////////////
	// total nr of points in each marginal
	int xres=GridToolsGetTotalPoints(muX.depth, muX.dimensions);
	int yres=GridToolsGetTotalPoints(muY.depth, muY.dimensions);

	// truncate marginals
	// both marginals are rounded to multiples of measureScale
	msg=MeasureToolsTruncateMeasures(muX.data, muY.data, xres, yres, measureScale);
	if(msg!=0) {
		return msg;
	}
	#endif

	
	
	///////////////////////////////////////////////
	// problem setup
	// as above	
	
	TMultiScaleSetupSingleGrid MultiScaleSetupX(&muX,depth);
	msg=MultiScaleSetupX.Setup();
	if(msg!=0) {
		printf("%d\n",msg);
		return msg;
	}


	TMultiScaleSetupSingleGrid MultiScaleSetupY(&muY,depth);
	msg=MultiScaleSetupY.Setup();
	if(msg!=0) {
		printf("%d\n",msg);
		return msg;
	}

	
	// for general W_p distances the shielding neighbourhoods can no longer be
	// constructed directly based on grid structure. extract more data for more general
	// shielding method
	MultiScaleSetupX.SetupGridNeighbours();
	MultiScaleSetupX.SetupRadii();
	MultiScaleSetupY.SetupRadii();
	
	// use different class than in W2Grid() to describe general p-th power of Euclidean distance
	TCostFunctionProvider_PEuclidean costFunctionProvider(
			MultiScaleSetupX.resH, MultiScaleSetupY.resH,
			MultiScaleSetupX.posH, MultiScaleSetupY.posH,
			MultiScaleSetupX.nLayers, MultiScaleSetupX.dim,
			p);

	
	// different, more involved ShieldGenerator class
	TShieldGeneratorTree_PEuclidean shieldGenerator(
			MultiScaleSetupX.dim,
			MultiScaleSetupY.HP, MultiScaleSetupY.posH, MultiScaleSetupY.radii,
			0 /* finest layer */, 0 /* coarsest layer */,
			MultiScaleSetupX.resH, MultiScaleSetupX.posH, MultiScaleSetupX.neighboursH,
			p,0.);
			
						
	TShortCutCouplingHandlerInterface couplingHandlerInterface(
			MultiScaleSetupX.resH, MultiScaleSetupY.resH,
			MultiScaleSetupX.nLayers);
	
	
	#ifdef USE_CPLEX
	// LP subsolver: CPLEX 
	TShortCutSubSolverInterfaceCPLEX subSolverInterface(
			MultiScaleSetupX.nLayers,
			MultiScaleSetupX.muH, MultiScaleSetupY.muH,
			&couplingHandlerInterface,true);
	
	static constexpr int VIOLATION_CHECKMODE=TShortCutSolver::VCHECK_DUAL;		
	
	#endif
	#ifdef USE_LEMON
	// LP subsolver: Lemon
	TShortCutSubSolverInterfaceLemon subSolverInterface(
			MultiScaleSetupX.nLayers,
			MultiScaleSetupX.muH, MultiScaleSetupY.muH,
			&couplingHandlerInterface,
			measureScale, cScale,
			dualOffset);

	static constexpr int VIOLATION_CHECKMODE=TShortCutSolver::VCHECK_PRIMAL;

	#endif
	#ifdef USE_LPSOLVE
	// LP subsolver: lp_solve 
	TShortCutSubSolverInterfaceLpSolve subSolverInterface(
			MultiScaleSetupX.nLayers,
			MultiScaleSetupX.muH, MultiScaleSetupY.muH,
			&couplingHandlerInterface,true);
	
	static constexpr int VIOLATION_CHECKMODE=TShortCutSolver::VCHECK_DUAL;		
	#endif
	
	
	///////////////////////////////////////////////

	TShortCutMultiScaleSolver MultiScaleSolver(
			&costFunctionProvider,
			&couplingHandlerInterface,
			&subSolverInterface,
			&shieldGenerator,
			MultiScaleSetupX.HP, MultiScaleSetupY.HP,
			1, // coarsest layer
			VIOLATION_CHECKMODE
			);
			
	MultiScaleSolver.autoDeletePointers=false;

	msg=MultiScaleSolver.solve();
	if(msg!=0) {
		printf("%d\n",msg);
		return msg;
	}
	
	printf("Wasserstein distance:%f\n",std::pow(MultiScaleSolver.objective,1/p));
	
	return 0;
	
}


int W2Grid_64x64() {
	int msg;
	
	char filenameMuX[]="data/grid/000_064.dat";
	char filenameMuY[]="data/grid/001_064.dat";

	std::vector<double> muXdat=readFile<double>(filenameMuX);
	std::vector<double> muYdat=readFile<double>(filenameMuY);

	
	// dimensions of grid
	int muXdim[]={64, 64};
	int muYdim[]={64, 64};
	
	
	// put raw marginal data into small container structs
	TDoubleMatrix muX,muY;

	muX.data=muXdat.data();
	muX.dimensions=muXdim;
	muX.depth=2;
	
	muY.data=muYdat.data();
	muY.dimensions=muYdim;
	muY.depth=2;
	
	// fundamental parameters
	int depth=6; // hierarchy depth

	#ifdef USE_LEMON
	// Lemon parameters	
	double measureScale=1E-9; // scale to which marginal measures are truncated
	double cScale=1E-7; // scale to which cost function is truncated
	bool dualOffset=true; // activate small trick with dual variables
		// that slightly accelerates the lemon network simplex
	
	///////////////////////////////////////////////
	// total nr of points in each marginal
	int xres=GridToolsGetTotalPoints(muX.depth, muX.dimensions);
	int yres=GridToolsGetTotalPoints(muY.depth, muY.dimensions);

	// truncate marginals
	// both marginals are rounded to multiples of measureScale
	msg=MeasureToolsTruncateMeasures(muX.data, muY.data, xres, yres, measureScale);
	if(msg!=0) {
		return msg;
	}
	#endif


	///////////////////////////////////////////////
	// problem setup
	
	
	// generate multi-scale problem representation:
	// the MultiScaleSetup object will subsequently hold all
	// required data for the solver and other algorithms
	
	// the class TMultiScaleSetupGrid assumes that muX and muY describe measures
	// that live on regular Cartesian grids with edge lenghts 1
	
	// the more general base class TMultiScaleSetupBasefor uses marginals with support on arbitrary point clouds
	TMultiScaleSetupSingleGrid MultiScaleSetupX(&muX,depth);
	msg=MultiScaleSetupX.Setup();
	if(msg!=0) {
		printf("%d\n",msg);
		return msg;
	}

	TMultiScaleSetupSingleGrid MultiScaleSetupY(&muY,depth);
	msg=MultiScaleSetupY.Setup();
	if(msg!=0) {
		printf("%d\n",msg);
		return msg;
	}

	
	// setup various aux components for actual solver algorithm
	TCostFunctionProvider_Dynamic costFunctionProvider(
			MultiScaleSetupX.resH, MultiScaleSetupY.resH,
			MultiScaleSetupX.posH, MultiScaleSetupY.posH,
			MultiScaleSetupX.nLayers, MultiScaleSetupX.dim);
	
	TShieldGeneratorGrid_SqrEuclidean shieldGenerator(
			MultiScaleSetupX.dim,
			MultiScaleSetupX.dimH, MultiScaleSetupY.dimH,
			MultiScaleSetupX.nLayers);
			
	TShortCutCouplingHandlerInterface couplingHandlerInterface(
			MultiScaleSetupX.resH, MultiScaleSetupY.resH,
			MultiScaleSetupX.nLayers);
	
	
	#ifdef USE_CPLEX
	// LP subsolver: CPLEX 
	TShortCutSubSolverInterfaceCPLEX subSolverInterface(
			MultiScaleSetupX.nLayers,
			MultiScaleSetupX.muH, MultiScaleSetupY.muH,
			&couplingHandlerInterface,true);
	
	static constexpr int VIOLATION_CHECKMODE=TShortCutSolver::VCHECK_DUAL;		
	
	#endif
	#ifdef USE_LEMON
	// LP subsolver: Lemon
	TShortCutSubSolverInterfaceLemon subSolverInterface(
			MultiScaleSetupX.nLayers,
			MultiScaleSetupX.muH, MultiScaleSetupY.muH,
			&couplingHandlerInterface,
			measureScale, cScale,
			dualOffset);

	static constexpr int VIOLATION_CHECKMODE=TShortCutSolver::VCHECK_PRIMAL;

	#endif
	#ifdef USE_LPSOLVE
	// LP subsolver: lp_solve 
	TShortCutSubSolverInterfaceLpSolve subSolverInterface(
			MultiScaleSetupX.nLayers,
			MultiScaleSetupX.muH, MultiScaleSetupY.muH,
			&couplingHandlerInterface,true);
	
	static constexpr int VIOLATION_CHECKMODE=TShortCutSolver::VCHECK_DUAL;		
	#endif


	///////////////////////////////////////////////
	// main solver object
	TShortCutMultiScaleSolver MultiScaleSolver(
			&costFunctionProvider,
			&couplingHandlerInterface,
			&subSolverInterface,
			&shieldGenerator,
			MultiScaleSetupX.HP, MultiScaleSetupY.HP,
			1, // coarsest layer
			VIOLATION_CHECKMODE
			);
			
	MultiScaleSolver.autoDeletePointers=false;

	msg=MultiScaleSolver.solve();
	if(msg!=0) {
		printf("%d\n",msg);
		return msg;
	}
	

	return 0;

}



int WpSphere() {
	int msg;
	
	// read raw data from file
	char filenamePoints[]="data/sphere/3081/points.dat";
	char filenameMuX[]="data/sphere/3081/m0.dat";
	char filenameMuY[]="data/sphere/3081/m1.dat";

	std::vector<double> posDat=readFile<double>(filenamePoints);
	std::vector<double> muXdat=readFile<double>(filenameMuX);
	std::vector<double> muYdat=readFile<double>(filenameMuY);

	
	// fundamental parameters
	int dim=3;
	int res=muXdat.size();
	int depth=5;
	int nNeighbours=5;
	double p=2;
	
	#ifdef USE_LEMON
	// Lemon parameters	
	double measureScale=1E-9; // scale to which marginal measures are truncated
	double cScale=1E-7; // scale to which cost function is truncated
	bool dualOffset=true; // activate small trick with dual variables
		// that slightly accelerates the lemon network simplex

	///////////////////////////////////////////////
	// truncate marginals
	// both marginals are rounded to multiples of measureScale
	msg=MeasureToolsTruncateMeasures(muXdat.data(), muYdat.data(), res, res, measureScale);
	#endif
	
	///////////////////////////////////////////////

	// setup problem data
	int posXdim[]={res, dim};
	
	TDoubleMatrix posX;

	posX.data=posDat.data();
	posX.dimensions=posXdim;
	posX.depth=2;
		
	
	
	///////////////////////////////////////////////

	TMultiScaleSetupSingleSphere MultiScaleSetupX(&posX,muXdat.data(),depth);
	MultiScaleSetupX.HierarchyBuilderChildMode=THierarchyBuilder::CM_Tree;
	msg=MultiScaleSetupX.Setup();	
	if(msg!=0) {
		printf("%d\n",msg);
		return msg;
	}


	TMultiScaleSetupSingleSphere MultiScaleSetupY(&posX,muYdat.data(),depth);
	MultiScaleSetupY.HierarchyBuilderChildMode=THierarchyBuilder::CM_Tree;
	msg=MultiScaleSetupY.Setup();	
	if(msg!=0) {
		printf("%d\n",msg);
		return msg;
	}

	
	printf("hierarchical cardinalities:\n");
	for(int layer=0;layer<MultiScaleSetupX.nLayers;layer++) {
		printf("%d\t%d\n",layer,MultiScaleSetupX.HP->layers[layer]->nCells);
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
	
	
	
	TCostFunctionProvider_Sphere costFunctionProvider(
			MultiScaleSetupX.resH, MultiScaleSetupY.resH,
			MultiScaleSetupX.posH, MultiScaleSetupY.posH,
			MultiScaleSetupX.nLayers, MultiScaleSetupX.dim,
			p);

	
			
	TShieldGeneratorTree_Sphere shieldGenerator(
			MultiScaleSetupX.dim,
			MultiScaleSetupY.HP, MultiScaleSetupY.posH, MultiScaleSetupY.radii,
			0 /* finest layer */, 0 /* coarsest layer */,
			MultiScaleSetupX.resH, MultiScaleSetupX.posH, MultiScaleSetupX.neighboursH,
			p);
	
	TShortCutCouplingHandlerInterface couplingHandlerInterface(
			MultiScaleSetupX.resH, MultiScaleSetupY.resH,
			MultiScaleSetupX.nLayers);
	
	#ifdef USE_CPLEX
	// LP subsolver: CPLEX 
	TShortCutSubSolverInterfaceCPLEX subSolverInterface(
			MultiScaleSetupX.nLayers,
			MultiScaleSetupX.muH, MultiScaleSetupY.muH,
			&couplingHandlerInterface,true);
	
	static constexpr int VIOLATION_CHECKMODE=TShortCutSolver::VCHECK_DUAL;		
	
	#endif
	#ifdef USE_LEMON
	// LP subsolver: Lemon
	TShortCutSubSolverInterfaceLemon subSolverInterface(
			MultiScaleSetupX.nLayers,
			MultiScaleSetupX.muH, MultiScaleSetupY.muH,
			&couplingHandlerInterface,
			measureScale, cScale,
			dualOffset);

	static constexpr int VIOLATION_CHECKMODE=TShortCutSolver::VCHECK_PRIMAL;

	#endif
	#ifdef USE_LPSOLVE
	// LP subsolver: lp_solve 
	TShortCutSubSolverInterfaceLpSolve subSolverInterface(
			MultiScaleSetupX.nLayers,
			MultiScaleSetupX.muH, MultiScaleSetupY.muH,
			&couplingHandlerInterface,true);
	
	static constexpr int VIOLATION_CHECKMODE=TShortCutSolver::VCHECK_DUAL;		
	#endif
	
	
	///////////////////////////////////////////////

	TShortCutMultiScaleSolver MultiScaleSolver(
			&costFunctionProvider,
			&couplingHandlerInterface,
			&subSolverInterface,
			&shieldGenerator,
			MultiScaleSetupX.HP, MultiScaleSetupY.HP,
			1, // coarsest layer
			VIOLATION_CHECKMODE
			);
			
	MultiScaleSolver.autoDeletePointers=false;

	msg=MultiScaleSolver.solve();
	if(msg!=0) {
		printf("%d\n",msg);
		return msg;
	}
	
	printf("Wasserstein distance:%f\n",std::pow(MultiScaleSolver.objective,1/p));
	
	
	return 0;

}

