#include"TCouplingHandler.h"

// Dense Data Structure

template<typename Vc>
TCouplingHandlerDensePrototype<Vc>::TCouplingHandlerDensePrototype(int _xres, int _yres, Vc *_c, double *_mu) {
	xres=_xres;
	yres=_yres;
	c=_c;
	mu=_mu;
	total=xres*yres;
}


template<typename Vc>
TCouplingHandlerDensePrototype<Vc>::~TCouplingHandlerDensePrototype() {
}


template<typename Vc>
void TCouplingHandlerDensePrototype<Vc>::clearMuRow(int x) {
	int y;
	for(y=0;y<yres;y++) {
		mu[yres*x+y]=0;
	}
}


template class TCouplingHandlerDensePrototype<double>;
template class TCouplingHandlerDensePrototype<int>;


// Semi-Dense Data Structure

// c and mu are stored in dense arrays, but only sparse lists of active variables are used.
// intention: can quickly alter 'active set' without significant reallocations.



template<typename Vc>
TCouplingHandlerSemiDensePrototype<Vc>::TCouplingHandlerSemiDensePrototype(int _xres, int _yres, Vc *_c, double *_mu, TVarListHandler *_xVars) {
	xres=_xres;
	yres=_yres;
	c=_c;
	mu=_mu;

	// assign xVars-pointer
	xVars=_xVars;
	// create yVars via transpose operation
	yVars= new TVarListHandler;
	yVars->fillViaTranspose(xVars,yres);
	total=xVars->total;
}


template<typename Vc>
TCouplingHandlerSemiDensePrototype<Vc>::~TCouplingHandlerSemiDensePrototype() {
	delete yVars;
}


template<typename Vc>
void TCouplingHandlerSemiDensePrototype<Vc>::clearMuRow(int x) {
	int y,yIndex;
	for(yIndex=0;yIndex<xVars->lenList->at(x);yIndex++) {
		y=xVars->varList[x]->at(yIndex);
		mu[yres*x+y]=0;
	}
}


template<typename Vc>
void TCouplingHandlerSemiDensePrototype<Vc>::updateXVars(TVarListHandler *newXVars, bool keepXVars) {
	// delete old x&y Vars
	if(!keepXVars) {
		delete xVars;
	}
	delete yVars;
	// assign xVars-pointer
	xVars=newXVars;
	total=xVars->total;
	// create yVars via transpose operation
	yVars= new TVarListHandler;
	yVars->fillViaTranspose(xVars,yres);
}


template class TCouplingHandlerSemiDensePrototype<double>;
template class TCouplingHandlerSemiDensePrototype<int>;


// Sparse Data Structure

// c and mu are stored in one contiguous block, just the data as in scipy.sparse.csr_matrix


TCouplingHandlerSparse::TCouplingHandlerSparse(int _xres, int _yres, TCostFunctionProviderBase *_cProvider,
		TVarListHandler *_xVars) {
	xres=_xres;
	yres=_yres;
	cProvider=_cProvider;

	// assign xVars-pointer
	xVars=_xVars;
	total=xVars->total;

	// initialize c
	c=cProvider->getC(xVars);
	free_c=cProvider->free_c();


	// allocate mu memory
	mu=(double*) malloc(sizeof(double)*total);
	
	// create offsets array
	offsets=(int*) malloc(sizeof(int)*xres);
	computeOffsets();
}

void TCouplingHandlerSparse::computeOffsets() {
	offsets[0]=0;
	for(int i=0;i<xres-1;i++) {
		offsets[i+1]=offsets[i]+xVars->lenList->at(i);
	}
}

TCouplingHandlerSparse::~TCouplingHandlerSparse() {
	free(mu);
	if(free_c) {
		free(c);
	}
	free(offsets);
}


void TCouplingHandlerSparse::clearMuRow(int x) {
	int yIndex;
	for(yIndex=0;yIndex<xVars->lenList->at(x);yIndex++) {
		mu[offsets[x]+yIndex]=0;
	}
}

void TCouplingHandlerSparse::updateXVars(TVarListHandler *newXVars, bool keepXVars) {

	double *muOld;
	TVarListHandler *xVarsOld;
	
	muOld=mu;
	xVarsOld=xVars;

	// assign xVars-pointer
	xVars=newXVars;
	total=xVars->total;
	
	// allocate mu memory
	mu=(double*) malloc(sizeof(double)*total);


	// copy old masses to new memory
	TVarListSignal<double> *signalMuOld = new TVarListSignal<double>(xVarsOld,muOld);
	TVarListSignal<double> *signalMuNew = new TVarListSignal<double>(xVars,mu);
	signalMuNew->transcribeSorted(signalMuOld,0.);
	delete signalMuOld;
	delete signalMuNew;

	// update cost array
	if(free_c) {
		free(c);
	}
	c=cProvider->getC(xVars);


	// delete old xVars
	if(!keepXVars) {
		delete xVarsOld;
	}
	// free old mu memory
	free(muOld);

	computeOffsets();
}


