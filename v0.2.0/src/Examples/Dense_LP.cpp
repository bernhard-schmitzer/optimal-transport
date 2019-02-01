#include<cstdlib>
#include <vector>

#include<Common.h>

#ifdef USE_CPLEX
#include<LP_CPLEX.h>
#endif
#ifdef USE_LEMON
#include<LP_Lemon.h>
#endif
#ifdef USE_LPSOLVE
#include<LP_lp_solve.h>
#endif

int main();
int W2Grid();



int main() {
	W2Grid();
}


int W2Grid() {
	int msg;
	
	// both marginals life on a regular 8x8 grid
	// provide measure data as consecutive 64 double entries
	double muXdat[]={6.101275e-03, 7.204110e-03, 8.535062e-03, 9.910280e-03, 1.104128e-02, 1.166015e-02, 1.167389e-02, 1.123175e-02, 6.704433e-03, 8.128397e-03, 9.854252e-03, 1.163492e-02, 1.309322e-02, 1.388739e-02, 1.390795e-02, 1.335771e-02, 7.292185e-03, 9.029484e-03, 1.114963e-02, 1.334184e-02, 1.513976e-02, 1.613150e-02, 1.620053e-02, 1.561280e-02, 7.908921e-03, 9.950044e-03, 1.245024e-02, 1.503625e-02, 1.715816e-02, 1.834562e-02, 1.848901e-02, 1.792387e-02, 8.567936e-03, 1.091527e-02, 1.379087e-02, 1.675992e-02, 1.919501e-02, 2.057699e-02, 2.080940e-02, 2.028992e-02, 9.247266e-03, 1.191367e-02, 1.517757e-02, 1.854434e-02, 2.131180e-02, 2.291030e-02, 2.324255e-02, 2.273922e-02, 9.914861e-03, 1.291900e-02, 1.659968e-02, 2.040437e-02, 2.355167e-02, 2.540296e-02, 2.582486e-02, 2.523923e-02, 1.055807e-02, 1.392742e-02, 1.806860e-02, 2.236741e-02, 2.594533e-02, 2.805923e-02, 2.849642e-02, 2.764292e-02};

	double muYdat[]={3.908070e-03, 4.687770e-03, 5.569132e-03, 6.429051e-03, 7.125903e-03, 7.561235e-03, 7.732784e-03, 7.742170e-03, 5.244827e-03, 5.949688e-03, 6.792079e-03, 7.684573e-03, 8.490201e-03, 9.083672e-03, 9.421825e-03, 9.576329e-03, 7.285414e-03, 7.776036e-03, 8.404373e-03, 9.162126e-03, 9.953290e-03, 1.064454e-02, 1.115002e-02, 1.149372e-02, 1.034929e-02, 1.044027e-02, 1.058721e-02, 1.093244e-02, 1.148517e-02, 1.214090e-02, 1.277317e-02, 1.333126e-02, 1.502886e-02, 1.456200e-02, 1.391201e-02, 1.346945e-02, 1.344750e-02, 1.383385e-02, 1.447324e-02, 1.520699e-02, 2.202719e-02, 2.099936e-02, 1.930211e-02, 1.767088e-02, 1.664565e-02, 1.640100e-02, 1.678885e-02, 1.752175e-02, 3.157802e-02, 3.026055e-02, 2.748538e-02, 2.437437e-02, 2.190360e-02, 2.055179e-02, 2.025589e-02, 2.062372e-02, 4.272636e-02, 4.166590e-02, 3.810900e-02, 3.349480e-02, 2.927364e-02, 2.633617e-02, 2.482542e-02, 2.433618e-02};
	
	// dimensions of grid
	int dim=2;
	int muXdim[]={8, 8};
	int muYdim[]={8, 8};
	
	
	// put raw marginal data into small container structs
	TDoubleMatrix muX,muY;

	muX.data=muXdat;
	muX.dimensions=muXdim;
	muX.depth=dim;
	
	muY.data=muYdat;
	muY.dimensions=muYdim;
	muY.depth=dim;
	
	
	TDoubleMatrix *posX=GridToolsGetGridMatrix(muX.depth, muX.dimensions);
	TDoubleMatrix *posY=GridToolsGetGridMatrix(muY.depth, muY.dimensions);
	
	int xres=posX->dimensions[0];
	int yres=posY->dimensions[0];
	std::vector<double> c(xres*yres);
	std::vector<double> mu(xres*yres);
	std::vector<double> alpha(xres);
	std::vector<double> beta(yres);
	
	#ifdef USE_LEMON
	// Lemon parameters	
	double measureScale=1E-9; // scale to which marginal measures are truncated
	double cScale=1E-7; // scale to which cost function is truncated
	
	// truncate marginals
	// both marginals are rounded to multiples of measureScale
	msg=MeasureToolsTruncateMeasures(muX.data, muY.data, xres, yres, measureScale);
	if(msg!=0) {
		return msg;
	}
	#endif

	
	
	TCostFunctionProvider_Dynamic costFunctionProvider(
		&xres, &yres,
		&(posX->data), &(posY->data),
		1, dim);

	costFunctionProvider.getCDense(c.data());
	
//	for(int i=0;i<yres;i++) {
//		printf(" %f",c[i]);
//	}
//	printf("\n\n");
	
	TCouplingHandlerDense couplingHandler(xres,yres,c.data(),mu.data());
	
	#ifdef USE_CPLEX
	TCPLEXNetSolver<TCouplingHandlerDense> solver(&couplingHandler,muX.data,muY.data,alpha.data(),beta.data());
	#endif
	#ifdef USE_LEMON
	TLemonSolver<TCouplingHandlerDense,TAlgorithm> solver(&couplingHandler,muX.data,muY.data,alpha.data(),beta.data(),
			measureScale,cScale);
	#endif
	#ifdef USE_LPSOLVE
	TLp_Solve<TCouplingHandlerDense> solver(&couplingHandler,muX.data,muY.data,alpha.data(),beta.data());
	#endif

	msg=solver.setup();
	if(msg!=0) {
		printf("error: %d\n",msg);
		return msg;
	}
	msg=solver.solve();
	if(msg!=0) {
		printf("error: %d\n",msg);
		return msg;
	}
	
	printf("objective: %e\n",solver.getObjective());

	
//	printf("optimal coupling:\n");
//	for(int x=0;x<xres;x++) {
//		for(int y=0;y<yres;y++) {
//			printf(" %f",mu[x*yres+y]);
//		}
//		printf("\n");
//	}

//	printf("dual variables:\n");
//	printf("alpha:\n");
//	for(int x=0;x<xres;x++) {
//		printf("\t%f\n",alpha[x]);
//	}
//	printf("beta:\n");
//	for(int y=0;y<yres;y++) {
//		printf("\t%f\n",beta[y]);
//	}


	
	freeTDoubleMatrix(posX);
	freeTDoubleMatrix(posY);

	
	return 0;

}




