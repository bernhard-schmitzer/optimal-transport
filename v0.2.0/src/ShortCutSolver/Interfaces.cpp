#include"Interfaces.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// TShortCutCouplingHandlerInterface
/////////////////////////////////////////////////////////////////////////////////////////////

TShortCutCouplingHandlerInterface::TShortCutCouplingHandlerInterface(int *_xresH, int *_yresH, int _nLayers) {
	xresH=_xresH;
	yresH=_yresH;
	nLayers=_nLayers;
	layer=0;
	couplingHandler=NULL;
	
}


TShortCutCouplingHandlerInterface::~TShortCutCouplingHandlerInterface() {
	if (couplingHandler!=NULL) {
		delete couplingHandler;
	}
}


int TShortCutCouplingHandlerInterface::setLayer(int newLayer,
		TCostFunctionProviderBase *costFunctionProvider, TVarListHandler *xVars) {

	if(newLayer>=nLayers) {
		return ERR_MULTISCALE_EXCEEDEDLEVELS;
	}
	if (couplingHandler!=NULL) {
		delete couplingHandler;
	}
	layer=newLayer;
	
	couplingHandler = new TCouplingHandlerSparse(xresH[layer],yresH[layer], costFunctionProvider, xVars);

	return 0;
}


int TShortCutCouplingHandlerInterface::getXres() {
	if (couplingHandler==NULL) { return 0; }
	return couplingHandler->xres;
	}

int TShortCutCouplingHandlerInterface::getYres() {
	if (couplingHandler==NULL) { return 0; }
	return couplingHandler->yres;
	}

TVarListHandler* TShortCutCouplingHandlerInterface::getXVars() {
	if (couplingHandler==NULL) { return NULL; }
	return couplingHandler->xVars;
	}



TVarListHandler* TShortCutCouplingHandlerInterface::getSupport() {
	if (couplingHandler==NULL) { return NULL; }

	int x,yIndex,y;
	TVarListHandler *result;
	result=new TVarListHandler();
	result->setupEmpty(couplingHandler->xres);
	for(x=0;x<couplingHandler->xres;x++) {
		for(yIndex=0;yIndex<couplingHandler->getRowLength(x);yIndex++) {
			if(couplingHandler->getMuRow(x,yIndex)>MASS_TOLERANCE) {
				y=couplingHandler->getColNr(x,yIndex);
				result->varList[x].push_back(y);
			}
		}
	}
	result->total=0;
	for(x=0;x<couplingHandler->xres;x++) {
		result->lenList[x]=result->varList[x].size();
		result->total+=result->lenList[x];
	}

	return result;
}

TVarListSignal<double>* TShortCutCouplingHandlerInterface::getSupportSignal() {
	if (couplingHandler==NULL) { return NULL; }

	int x,yIndex,y;
	TVarListHandler *varList;
	std::vector<double> signal(0);
	double *signalArray;
	TVarListSignal<double> *result;

	double mass;

	varList=new TVarListHandler();
	varList->setupEmpty(couplingHandler->xres);

	for(x=0;x<couplingHandler->xres;x++) {
		for(yIndex=0;yIndex<couplingHandler->getRowLength(x);yIndex++) {
			mass=couplingHandler->getMuRow(x,yIndex);
			if(mass>MASS_TOLERANCE) {
				y=couplingHandler->getColNr(x,yIndex);
				varList->varList[x].push_back(y);
				signal.push_back(mass);
			}
		}
	}

	// setup varList counting variables
	varList->total=0;
	for(x=0;x<couplingHandler->xres;x++) {
		varList->lenList[x]=varList->varList[x].size();
		varList->total+=varList->lenList[x];
	}

	// setup signal array
	signalArray=(double*) malloc(sizeof(double)*signal.size());
	std::copy(signal.begin(), signal.end(), signalArray);

	result=new TVarListSignal<double>(varList,signalArray);

	return result;
}


TSparseCSRContainer TShortCutCouplingHandlerInterface::getCouplingCSR() {
	TSparseCSRContainer result;
	if (couplingHandler==NULL) { return result; }

	int x,yIndex,y;
	int count=0;

	double mass;

	result.xres=couplingHandler->xres;
	result.yres=couplingHandler->yres;
	result.indptr.resize(result.xres+1);
	result.indptr[0]=0;

	for(x=0;x<couplingHandler->xres;x++) {
		for(yIndex=0;yIndex<couplingHandler->getRowLength(x);yIndex++) {
			mass=couplingHandler->getMuRow(x,yIndex);
			if(mass>MASS_TOLERANCE) {
				y=couplingHandler->getColNr(x,yIndex);
				result.indices.push_back(y);
				result.data.push_back(mass);
				count++;
			}
		}
		result.indptr[x+1]=count;
	}

	result.nonZeros=count;
	return result;
}


TSparsePosContainer TShortCutCouplingHandlerInterface::getCouplingPos() {
	TSparsePosContainer result;
	if (couplingHandler==NULL) { return result; }

	int x,yIndex,y;
	int count=0;

	double mass;

	result.xres=couplingHandler->xres;
	result.yres=couplingHandler->yres;

	for(x=0;x<couplingHandler->xres;x++) {
		for(yIndex=0;yIndex<couplingHandler->getRowLength(x);yIndex++) {
			mass=couplingHandler->getMuRow(x,yIndex);
			if(mass>MASS_TOLERANCE) {
				y=couplingHandler->getColNr(x,yIndex);
				result.posStart.push_back(x);
				result.posEnd.push_back(y);
				result.mass.push_back(mass);
				count++;
			}
		}
	}

	result.nParticles=count;
	return result;
}


void TShortCutCouplingHandlerInterface::updateXVars(TVarListHandler *newXVars, bool keepXVars) {
	if (couplingHandler==NULL) { return; }

	couplingHandler->updateXVars(newXVars, keepXVars);
}


bool TShortCutCouplingHandlerInterface::dualViolationCheck(bool doProjection, double *alpha, double *beta) {
	if (couplingHandler==NULL) { return true; }

	// check dual variables (available through solverInterface->alpha and ->beta)
	// for violation of relevant dual constraints (list of those available via couplingHandler)
	// if doProjection=true: reduce alpha to match constraints and clear corresponding muRow
	//     this is useful, for example, to re-initialize the hungarian method
	// if doProjection=false: after first found violated constraint just return true.
	//     useful, for CPLEX, where dual variables are just a means of checking optimality
	bool result;
	double alphaMax;
	int x,yIndex,y;
	result=false;
	for(x=0;x<couplingHandler->xres;x++) {
		alphaMax=SLACK_INFINITY;
		for(yIndex=0;yIndex<couplingHandler->getRowLength(x);yIndex++) {
			y=couplingHandler->getColNr(x,yIndex);
			alphaMax=std::min(alphaMax,couplingHandler->getCRow(x,yIndex)-beta[y]);
		}
		if(alpha[x]>alphaMax+SLACK_TOLERANCE) {
				if(!doProjection) {
					return true;
				}
				alpha[x]=alphaMax;
				couplingHandler->clearMuRow(x);
				result=true;
		}
	}
	return result;
}



/////////////////////////////////////////////////////////////////////////////////////////////
// TSubSolverInterface dummy methods
/////////////////////////////////////////////////////////////////////////////////////////////

int TShortCutSubSolverInterfaceBase::solve() {return ERR_BASE_NOTIMPLEMENTED;}
int TShortCutSubSolverInterfaceBase::prepareUpdate(__attribute__((unused)) TVarListHandler *newXVars) {
	return ERR_BASE_NOTIMPLEMENTED;
}
int TShortCutSubSolverInterfaceBase::update() {return ERR_BASE_NOTIMPLEMENTED;}
double TShortCutSubSolverInterfaceBase::getObjective() {return 0.;}



