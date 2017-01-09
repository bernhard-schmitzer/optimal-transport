#include<stdlib.h>
#include<PythonTypes.h>

#include"TCouplingHandler.h"
#include"TCPLEXNetSolver.h"

using namespace std;

template class TCPLEXNetSolver<TCouplingHandlerDense>;
template class TCPLEXNetSolver<TCouplingHandlerSemiDense>;
template class TCPLEXNetSolver<TCouplingHandlerSparse>;


extern "C" {

int solveDense(TDoubleMatrix *c, TDoubleMatrix *muX, TDoubleMatrix *muY,
	TDoubleMatrix *mu, bool getDualVariables, TDoubleMatrix *alpha, TDoubleMatrix *beta) {

	int xres,yres,msg;

	xres=c->dimensions[0];
	yres=c->dimensions[1];

	TCouplingHandlerDense couplingHandler(xres,yres,c->data,mu->data);
	TCPLEXNetSolver<TCouplingHandlerDense> *solver;

	if (getDualVariables) {
		solver = new TCPLEXNetSolver<TCouplingHandlerDense>(&couplingHandler, muX->data, muY->data,
				alpha->data, beta->data);
	} else {
		solver = new TCPLEXNetSolver<TCouplingHandlerDense>(&couplingHandler, muX->data, muY->data);
	}

	msg=solver->setup();
	if(msg!=0) {
		return msg;
	}

	msg=solver->solve();

	delete solver;
	return msg;
}

int solveSemiDense(TDoubleMatrix *c, TDoubleMatrix *muX, TDoubleMatrix *muY,
	TDoubleMatrix *mu,
	bool getDualVariables, TDoubleMatrix *alpha, TDoubleMatrix *beta,
	TInteger32Matrix *VLXindices, TInteger32Matrix *VLXindptr
	) {

	int xres,yres,msg;

	xres=c->dimensions[0];
	yres=c->dimensions[1];

	TVarListHandler *xVars;
	xVars=new TVarListHandler();
	xVars->fillFromCSRIndexList(VLXindices->data, VLXindptr->data, xres, VLXindices->dimensions[0]);

	TCouplingHandlerSemiDense couplingHandler(xres,yres,c->data,mu->data,xVars);

	TCPLEXNetSolver<TCouplingHandlerSemiDense> *solver;

	if (getDualVariables) {
		solver = new TCPLEXNetSolver<TCouplingHandlerSemiDense>(&couplingHandler, muX->data, muY->data,
				alpha->data,beta->data);
	} else {
		solver = new TCPLEXNetSolver<TCouplingHandlerSemiDense>(&couplingHandler, muX->data, muY->data);
	}

	msg=solver->setup();
	if(msg!=0) {
		return msg;
	}

	msg=solver->solve();

	// Clean up memory
	delete xVars;
	delete solver;

	return msg;
}

int solveSparse(int xres, int yres, TDoubleMatrix *c, TDoubleMatrix *muX, TDoubleMatrix *muY,
	bool getDualVariables, TDoubleMatrix *alpha, TDoubleMatrix *beta,
	TInteger32Matrix *VLXindices, TInteger32Matrix *VLXindptr,
	TInteger32Matrix *muSpecs, TInteger64Matrix *couplingHandlerAddr
	) {

	int msg;

	TVarListHandler *xVars;
	xVars=new TVarListHandler();
	xVars->fillFromCSRIndexList(VLXindices->data, VLXindptr->data, xres, VLXindices->dimensions[0]);

	TCostFunctionProviderStatic cProvider(c->data);
	TCouplingHandlerSparse *couplingHandler = new TCouplingHandlerSparse(xres,yres, &cProvider, xVars);

	TCPLEXNetSolver<TCouplingHandlerSparse> *solver;

	if (getDualVariables) {
		solver = new TCPLEXNetSolver<TCouplingHandlerSparse>(couplingHandler, muX->data, muY->data,
				alpha->data,beta->data);
	} else {
		solver = new TCPLEXNetSolver<TCouplingHandlerSparse>(couplingHandler, muX->data, muY->data);
	}

	msg=solver->setup();
	if(msg!=0) {
		return msg;
	}

	msg=solver->solve();

	muSpecs->data[0]=couplingHandler->total;
	couplingHandlerAddr->data[0]=(long int) couplingHandler;
	// Clean up memory
	//delete xVars;
	delete solver;

	return msg;
}

int solveSparse_fullC(TDoubleMatrix *c, TDoubleMatrix *muX, TDoubleMatrix *muY,
	bool getDualVariables, TDoubleMatrix *alpha, TDoubleMatrix *beta,
	TInteger32Matrix *VLXindices, TInteger32Matrix *VLXindptr,
	TInteger32Matrix *muSpecs, TInteger64Matrix *couplingHandlerAddr
	) {

	int xres,yres,msg;

	xres=c->dimensions[0];
	yres=c->dimensions[1];

	TVarListHandler *xVars;
	xVars=new TVarListHandler();
	xVars->fillFromCSRIndexList(VLXindices->data, VLXindptr->data, xres, VLXindices->dimensions[0]);

	TCostFunctionProviderFullArray cProvider(xres,yres,c->data);
	TCouplingHandlerSparse *couplingHandler = new TCouplingHandlerSparse(xres,yres, &cProvider, xVars);

	TCPLEXNetSolver<TCouplingHandlerSparse> *solver;

	if (getDualVariables) {
		solver = new TCPLEXNetSolver<TCouplingHandlerSparse>(couplingHandler, muX->data, muY->data,
				alpha->data,beta->data);
	} else {
		solver = new TCPLEXNetSolver<TCouplingHandlerSparse>(couplingHandler, muX->data, muY->data);
	}

	msg=solver->setup();
	if(msg!=0) {
		return msg;
	}

	msg=solver->solve();

	muSpecs->data[0]=couplingHandler->total;
	couplingHandlerAddr->data[0]=(long int) couplingHandler;
	// Clean up memory
	delete solver;

	return msg;
}

int solveSparse_getMu(TDoubleMatrix *mu, long int couplingHandlerAddr) {
	TCouplingHandlerSparse *couplingHandler = (TCouplingHandlerSparse*) couplingHandlerAddr;
	for(int i=0;i<mu->dimensions[0];i++) {
			mu->data[i]=couplingHandler->mu[i];
	}
	delete couplingHandler->xVars;
	delete couplingHandler;
	return 0;
}

}

