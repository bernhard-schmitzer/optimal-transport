#ifndef TCPLEXNetSolver_H_
#define TCPLEXNetSolver_H_

#include <cstdlib>
#include <ilcplex/cplexx.h>
#include <Common/ErrorCodes.h>
#include <Common/TCouplingHandler.h>
#include <Common/Tools.h>


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

	CPXENVptr CPXenv;
	CPXNETptr CPXnet;
	bool setupStatus; // true if environment is setup, false else
	bool solutionStatus;
	bool useDualVariables;
	bool muYFlippedStatus;

	double *muX, *muY;
	double *muYFlipped; // negative scaled copy of muY to model mass sinks
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
	
	void setupMuYFlipped();


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
	setupMuYFlipped(); // setup negatively scaled y marginal
}

template <class TCouplingHandlerType>
TCPLEXNetSolver<TCouplingHandlerType>::TCPLEXNetSolver(TCouplingHandlerType *_CouplingHandler,
		double *_muX, double *_muY, double *_alpha, double *_beta) :
		TCPLEXNetSolverBase(_muX, _muY, _alpha, _beta) {
	CouplingHandler=_CouplingHandler;
	xres=CouplingHandler->xres;
	yres=CouplingHandler->yres;
	setupMuYFlipped(); // setup negatively scaled y marginal
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
		// set xside, yside and cost function for all arcs
		for(yindex=0;yindex<rowlen;yindex++) {
			xrange[yindex]=xpos;
			yrange[yindex]=xres+CouplingHandler->getColNr(xpos,yindex); // do not forget offset by xres nr of nodes
			crange[yindex]=CouplingHandler->getCRow(xpos,yindex);
		}
		CPXXNETaddarcs(CPXenv, CPXnet, rowlen, xrange, yrange, NULL, NULL, crange, NULL);
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

	status=CPXXNETprimopt(CPXenv, CPXnet);
	if(status!=0) {
		return RSLT_ERROR_NOOPT;
	}

	// allocate space for optimal coupling
	double *mu;
	mu=(double*) malloc(sizeof(double)*nArcs);

	CPXXNETsolution(CPXenv, CPXnet, &status, &objective, mu, NULL, NULL, NULL);

	// checking status, sending abort in case something went wrong.
	if(status!=CPX_STAT_OPTIMAL) {
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
	
	
	// do marginal check
	//double marginalError=((TCouplingHandlerSparse*) CouplingHandler)->checkMarginalConstraints(muX,muY);
	//eprintf("marg error: %e\n",marginalError);

	// transfer optimal dual variables, if requested
	if(useDualVariables) {
		extractDuals();
	}

	solutionStatus=true;

	return 0;

}


#endif
