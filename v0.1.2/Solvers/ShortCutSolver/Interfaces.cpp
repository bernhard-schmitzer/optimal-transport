#include"Interfaces.h"

/////////////////////////////////////////////////////////////////////////////////////////////
// Dummy TCPLEXNetSolverInterface Methods
/////////////////////////////////////////////////////////////////////////////////////////////

int TSolverInterface::solve() {return MSG_NOT_IMPLEMENTED;}
int TSolverInterface::prepareUpdate() {return MSG_NOT_IMPLEMENTED;}
int TSolverInterface::update() {return MSG_NOT_IMPLEMENTED;}
double TSolverInterface::getObjective() {return 0.;}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TCouplingHandlerExtBase::TCouplingHandlerExtBase() {}

TCouplingHandlerExtBase::~TCouplingHandlerExtBase() {}


int TCouplingHandlerExtBase::getXres() {return 0;}
int TCouplingHandlerExtBase::getYres() {return 0;}
TVarListHandler* TCouplingHandlerExtBase::getXVars() { return NULL;}
TVarListHandler* TCouplingHandlerExtBase::getSupport() {return NULL;}
TVarListSignal<double>* TCouplingHandlerExtBase::getSupportSignal() {return NULL;}
void TCouplingHandlerExtBase::updateXVars(__attribute__((unused)) TVarListHandler *_newXVars,
		__attribute__((unused)) bool keepXVars) {}
bool TCouplingHandlerExtBase::dualViolationCheck(__attribute__((unused)) bool doProjection,
		__attribute__((unused)) double *alpha, __attribute__((unused)) double *beta) {return false;}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class TCouplingHandlerType>
TCouplingHandlerExt<TCouplingHandlerType>::TCouplingHandlerExt() {
	couplingHandler=NULL;
	deleteHandlerOnDestroy=false;
}


template<class TCouplingHandlerType>
TCouplingHandlerExt<TCouplingHandlerType>::TCouplingHandlerExt(TCouplingHandlerType *_couplingHandler,
		bool _deleteHandlerOnDestroy) {
	couplingHandler=_couplingHandler;
	deleteHandlerOnDestroy=_deleteHandlerOnDestroy;
}

template<class TCouplingHandlerType>
TCouplingHandlerExt<TCouplingHandlerType>::~TCouplingHandlerExt() {
	if (deleteHandlerOnDestroy && (couplingHandler!=NULL)) {
		delete couplingHandler;
	}
}


template<class TCouplingHandlerType>
int TCouplingHandlerExt<TCouplingHandlerType>::getXres() { return couplingHandler->xres; }

template<class TCouplingHandlerType>
int TCouplingHandlerExt<TCouplingHandlerType>::getYres() { return couplingHandler->yres; }

template<class TCouplingHandlerType>
TVarListHandler* TCouplingHandlerExt<TCouplingHandlerType>::getXVars() { return couplingHandler->xVars; }



template<class TCouplingHandlerType>
TVarListHandler* TCouplingHandlerExt<TCouplingHandlerType>::getSupport() {
	int x,yIndex,y;
	TVarListHandler *result;
	result=new TVarListHandler();
	result->setupEmpty(couplingHandler->xres);
	for(x=0;x<couplingHandler->xres;x++) {
		//cout << x << " " << (int) (couplingHandler->getRowLength(x)) << endl;
		for(yIndex=0;yIndex<couplingHandler->getRowLength(x);yIndex++) {
			if(couplingHandler->getMuRow(x,yIndex)>MASS_TOLERANCE) {
				y=couplingHandler->getColNr(x,yIndex);
				result->varList[x]->push_back(y);
			}
		}
	}
	result->total=0;
	for(x=0;x<couplingHandler->xres;x++) {
		result->lenList->at(x)=result->varList[x]->size();
		result->total+=result->lenList->at(x);
	}

	return result;
}

template<class TCouplingHandlerType>
TVarListSignal<double>* TCouplingHandlerExt<TCouplingHandlerType>::getSupportSignal() {
	int x,yIndex,y;
	TVarListHandler *varList;
	vector<double> signal(0);
	double *signalArray;
	TVarListSignal<double> *result;

	double mass;

	varList=new TVarListHandler();
	varList->setupEmpty(couplingHandler->xres);

	for(x=0;x<couplingHandler->xres;x++) {
		//cout << x << " " << (int) (couplingHandler->getRowLength(x)) << endl;
		for(yIndex=0;yIndex<couplingHandler->getRowLength(x);yIndex++) {
			mass=couplingHandler->getMuRow(x,yIndex);
			if(mass>MASS_TOLERANCE) {
				y=couplingHandler->getColNr(x,yIndex);
				varList->varList[x]->push_back(y);
				signal.push_back(mass);
			}
		}
	}

	// setup varList counting variables
	varList->total=0;
	for(x=0;x<couplingHandler->xres;x++) {
		varList->lenList->at(x)=varList->varList[x]->size();
		varList->total+=varList->lenList->at(x);
	}

	// setup signal array
	signalArray=(double*) malloc(sizeof(double)*signal.size());
	copy(signal.begin(), signal.end(), signalArray);

	result=new TVarListSignal<double>(varList,signalArray);

	return result;
}



template<class TCouplingHandlerType>
void TCouplingHandlerExt<TCouplingHandlerType>::updateXVars(TVarListHandler *newXVars, bool keepXVars) {
	couplingHandler->updateXVars(newXVars, keepXVars);
}


template<class TCouplingHandlerType>
bool TCouplingHandlerExt<TCouplingHandlerType>::dualViolationCheck(bool doProjection, double *alpha, double *beta) {
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
			alphaMax=min(alphaMax,couplingHandler->getCRow(x,yIndex)-beta[y]);
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


template class TCouplingHandlerExt<TCouplingHandlerSemiDense>;
template class TCouplingHandlerExt<TCouplingHandlerSparse>;
