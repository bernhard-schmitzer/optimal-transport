#include"TEpsScaling.h"

TEpsScalingHandler::TEpsScalingHandler(double _epsStart, double _epsTarget, int _epsSteps) {
	epsStart=_epsStart;
	epsTarget=_epsTarget;
	epsSteps=_epsSteps;
	
	epsList=(double*) malloc(sizeof(double)*(epsSteps+1));
	
	epsList[0]=epsStart;
	epsList[epsSteps]=epsTarget;
	
	for(int i=1;i<epsSteps;i++) {
		epsList[i]=pow(epsStart,1.-(double)i/epsSteps)*pow(epsTarget,(double)i/epsSteps);
	}
	
	epsScales=NULL;
	nLayers=0;
	epsLists=NULL;
	nEpsLists=NULL;
}

TEpsScalingHandler::~TEpsScalingHandler() {
	free(epsList);
	if(epsScales!=NULL) {
		free(epsScales);
	}
	if(nEpsLists!=NULL) {
		for(int i=0;i<nLayers;i++) {
			if(epsLists[i]!=NULL) {
				free(epsLists[i]);
			}
		}
		free(epsLists);
		free(nEpsLists);
	}
}

void TEpsScalingHandler::getEpsScalesFromBox(double boxScale, double layerExponent, int _nLayers) {
	nLayers=_nLayers;
	epsScales=(double*) malloc(sizeof(double)*nLayers);
	for(int i=0;i<nLayers-1;i++) {
		epsScales[i]=pow(boxScale*pow(0.5,i),layerExponent);
	}
	epsScales[nLayers-1]=0;
}


int TEpsScalingHandler::getEpsScalingSplit(int levelCoarsest, int nOverlap) {


	std::vector<std::vector<double> > vecEpsLists(nLayers);

	// divide eps list over scales
	int currentScale=levelCoarsest;
	for(int i=0;i<(epsSteps+1);i++) {
		double eps=epsList[i];
		while(eps<epsScales[currentScale]) {
			currentScale++;
		}
		vecEpsLists[currentScale].push_back(eps);
	}

	// create overlaps
	for(int i=1;i<nLayers;i++) {
		if(vecEpsLists[i-1].size()>0) {
			// if prev eps list is not empty
			int nPrev=vecEpsLists[i-1].size();
			int nOl=std::min(nPrev,nOverlap); // add at most nOverlap or number of elements
			for(int j=0;j<nOl;j++) {
				vecEpsLists[i].insert(vecEpsLists[i].begin(), vecEpsLists[i-1][nPrev-nOl+j]);
			}
		}
	}

	// check if coarsest eps list has at least one element
	if(vecEpsLists[levelCoarsest].size()==0) {
		return ERR_EPSSCALING_EMPTYCOARSESUBLIST;
	}
	
	// convert eps list to oldschool array format
	nEpsLists=(int*) malloc(sizeof(int)*nLayers);
	epsLists=(double**) malloc(sizeof(double*)*nLayers);
	for(int i=0;i<nLayers;i++) {
		nEpsLists[i]=vecEpsLists[i].size();
		if(nEpsLists[i]>0) {
			epsLists[i]=(double*) malloc(sizeof(double)*(nEpsLists[i]));
			for(int j=0;j<nEpsLists[i];j++) {
				epsLists[i][j]=vecEpsLists[i][j];
			}
		} else {
			epsLists[i]=NULL;
		}
	}
	
	return 0;
	
}


void TEpsScalingHandler::getEpsScalingAllFinest(int levelFinest) {
	nLayers=levelFinest+1;
	
	
	nEpsLists=(int*) malloc(sizeof(int)*nLayers);
	epsLists=(double**) malloc(sizeof(double*)*nLayers);
	
	// initialize all eps lists above finest layer as empty
	for(int i=0;i<nLayers-1;i++) {
		nEpsLists[i]=0;
		epsLists[i]=NULL;
	}
	
	// set finest layer list to standard list
	nEpsLists[nLayers-1]=epsSteps+1;
	epsLists[nLayers-1]=(double*) malloc(sizeof(double)*(epsSteps+1));
	memcpy(epsLists[nLayers-1],epsList,sizeof(double)*(epsSteps+1));
	
}

