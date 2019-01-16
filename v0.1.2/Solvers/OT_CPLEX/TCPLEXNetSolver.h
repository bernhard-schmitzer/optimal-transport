#ifndef TCPLEXNetSolver_H_
#define TCPLEXNetSolver_H_

#include<stdlib.h>
#include <ilcplex/cplex.h>
#include"TCouplingHandler.h"
#include<tools.h>


using namespace std;

class TCPLEXNetSolverBase {
public:
	static const int RSLT_ERROR_NOOPT=-101;
	static const int RSLT_ERROR_NOSOL=-102;
	static const int RSLT_ERROR_FREEPROB=-103;
	static const int RSLT_ERROR_FREEENV=-104;

	static const int SETUP_ERROR_NOSETUP=-201;
	static const int SETUP_ERROR_ALREADYSET=-202;
	static const int SETUP_ERROR_ENV=-203;
	static const int SETUP_ERROR_NET=-204;

	static const int SETUP_ERROR_NODUALS=-205;
	static const int SETUP_ERROR_EXTRACTDUALS=-206;

	static const int MSG_NOT_IMPLEMENTED = -207;

	CPXENVptr CPXenv;
	CPXNETptr CPXnet;
	bool setupStatus; // true if environment is setup, false else
	bool solutionStatus;
	bool useDualVariables;

	double *muX, *muY;
	double *alpha, *beta;
	double objective;
	int xres, yres, nArcs;

	TCPLEXNetSolverBase(double *_muX, double *_muY);
	TCPLEXNetSolverBase(double *_muX, double *_muY, double *_alpha, double *_beta);

	virtual ~TCPLEXNetSolverBase();
	virtual int setupEnvironment();
	virtual int setupMarginals();
	virtual int setupArcs();
	virtual int setupBase(int *arcStatus, int *nodeStatus);

	virtual int setup();

	virtual int solve();

	virtual int extractDuals();
	virtual int extractBase(int *arcStatus, int *nodeStatus);

	virtual int deleteArcs();

	virtual int cleanup();
	virtual double getObjective();


};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <class TCouplingHandlerType>
class TCPLEXNetSolver : public TCPLEXNetSolverBase {
public:

	TCouplingHandlerType *CouplingHandler;


	TCPLEXNetSolver(TCouplingHandlerType *_CouplingHandler, double *_muX, double *_muY);
	TCPLEXNetSolver(TCouplingHandlerType *_CouplingHandler, double *_muX, double *_muY, double *_alpha, double *_beta);

	int setupArcs();

	int solve();


};






template <class TCouplingHandlerType>
TCPLEXNetSolver<TCouplingHandlerType>::TCPLEXNetSolver(TCouplingHandlerType *_CouplingHandler,
		double *_muX, double *_muY) : TCPLEXNetSolverBase(_muX, _muY) {
	CouplingHandler=_CouplingHandler;
	xres=CouplingHandler->xres;
	yres=CouplingHandler->yres;
}

template <class TCouplingHandlerType>
TCPLEXNetSolver<TCouplingHandlerType>::TCPLEXNetSolver(TCouplingHandlerType *_CouplingHandler,
		double *_muX, double *_muY, double *_alpha, double *_beta) :
		TCPLEXNetSolverBase(_muX, _muY, _alpha, _beta) {
	CouplingHandler=_CouplingHandler;
	xres=CouplingHandler->xres;
	yres=CouplingHandler->yres;
}

template <class TCouplingHandlerType>
int TCPLEXNetSolver<TCouplingHandlerType>::setupArcs() {
	if (!setupStatus) {
		return SETUP_ERROR_NOSETUP;
	}

	int xpos,yindex,rowlen;
	int *xrange, *yrange;
	double *crange;

	xrange = (int*) malloc(sizeof(int)*yres);
	yrange = (int*) malloc(sizeof(int)*yres);
	crange = (double*) malloc(sizeof(double)*yres);

	// iterate through all rows
	for(xpos=0;xpos<xres;xpos++) {
		// for each row, iterate through columns of that row
		rowlen=CouplingHandler->getRowLength(xpos);
		//cout << "\tx: " << xpos << " rowlen: " << rowlen << endl;
		// set xside, yside and cost function for all arcs
		for(yindex=0;yindex<rowlen;yindex++) {
			//cout << "\t\t" << yindex << endl;
			xrange[yindex]=xpos;
			yrange[yindex]=xres+CouplingHandler->getColNr(xpos,yindex); // do not forget offset by xres nr of nodes
			//cout << "\t\t" << yrange[yindex] << endl;
			crange[yindex]=CouplingHandler->getCRow(xpos,yindex);
			//cout << "\t\t" << crange[yindex] << endl;
		}
		CPXNETaddarcs(CPXenv, CPXnet, rowlen, xrange, yrange, NULL, NULL, crange, NULL);
		nArcs+=(int)rowlen;
	}

	free(xrange);
	free(yrange);
	free(crange);

	solutionStatus=false;
	return 0;
}



template <class TCouplingHandlerType>
int TCPLEXNetSolver<TCouplingHandlerType>::solve() {
	int status;

	status=CPXNETprimopt(CPXenv, CPXnet);
	if(status!=0) {
		return RSLT_ERROR_NOOPT;
	}

	// allocate space for optimal coupling
	double *mu;
	mu=(double*) malloc(sizeof(double)*nArcs);

	CPXNETsolution(CPXenv, CPXnet, &status, &objective, mu, NULL, NULL, NULL);

	// checking status, sending abort in case something went wrong.
	if(status!=CPX_STAT_OPTIMAL) {
		//cout << "error orccured during solution retrieval: " << status << endl;
		free(mu);
		return RSLT_ERROR_NOSOL;
	}

	// transfer optimal coupling to CouplingHandlerObject
	int xpos,offset;
	offset=0;
	for(xpos=0;xpos<xres;xpos++) {
		CouplingHandler->setMuRow(xpos,mu+offset);
		offset+=CouplingHandler->getRowLength(xpos);
	}
	// cleaning up
	free(mu);

	// transfer optimal dual variables, if requested
	if(useDualVariables) {
		extractDuals();
	}

	solutionStatus=true;

	return 0;

}


#endif
