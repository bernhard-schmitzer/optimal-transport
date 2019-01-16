#include "TVarListHandler.h"

TVarListHandler::TVarListHandler() {
	res=0;
	total=0;
	lenList=NULL;
	varList=NULL;
}

// copy constructor
TVarListHandler::TVarListHandler(TVarListHandler *base) {
	res=base->res;
	total=base->total;
	lenList=new vector<int>(*(base->lenList));
	varList=NULL;

	varList=(vector<int>**) malloc(sizeof(vector<int>*)*res);
	for(int x=0;x<res;x++) {
		varList[x]=new vector<int>(*(base->varList[x]));
	}
}



TVarListHandler::~TVarListHandler() {
	clear();
}

void TVarListHandler::clear() {
	if(lenList!=NULL) {
			for(int i=0;i<res;i++) {
				delete varList[i];
			}
			free(varList);
			delete lenList;
	}
	varList=NULL;
	lenList=NULL;
	res=0;
	total=0;
}

void TVarListHandler::setupEmpty(int _res) {
	res=_res;
	lenList=new vector<int>(res);
	varList=(vector<int>**) malloc(sizeof(vector<int>*)*res);
	for(int x=0;x<res;x++) {
		varList[x]=new vector<int>(0);
	}
}

void TVarListHandler::fillViaTranspose(TVarListHandler *transpose, int yres) {
	int x,y,yIndex;
	clear();
	setupEmpty(yres);

	for(x=0;x<transpose->res;x++) {
		for(yIndex=0;yIndex<transpose->lenList->at(x);yIndex++) {
			y=transpose->varList[x]->at(yIndex);
			varList[y]->push_back(x);
		}
	}

	total=0;
	for(y=0;y<yres;y++) {
		lenList->at(y)=varList[y]->size();
		total+=lenList->at(y);
	}

}

void TVarListHandler::fillFromCSRIndexList(int *indices, int *indptr, int _res, int _total) {
	setupEmpty(_res);
	total=_total;
	int rowLen,offset;
	int x,y;
	for(x=0;x<_res;x++) {
		rowLen=indptr[x+1]-indptr[x];
		offset=indptr[x];
		(*lenList)[x]=(int) rowLen;
		varList[x]->resize(rowLen);
		for(y=0;y<rowLen;y++) {
			(*(varList[x]))[y]=indices[offset+y];
		}
	}
}

void TVarListHandler::writeToCSRIndexList(int *indices, int *indptr) {
	// write varList data to CSR (as in scipy.sparse.csr_matrix) list. indices must have size total, indptr must have size res+1
	int x,yIndex,offset;
	indptr[0]=0;
	offset=0;
	for(x=0;x<res;x++) {
		for(yIndex=0;yIndex<(*lenList)[x];yIndex++) {
			indices[offset]=(*(varList[x]))[yIndex];
			offset++;
		}
		indptr[x+1]=offset;
	}
	
}


void TVarListHandler::addToLine(int x, vector<int>* yCandidates) {
	int yIndex1,yIndex2,y;
	bool duplicate;
	for(yIndex1=0;yIndex1<(int) yCandidates->size();yIndex1++) {
		y=yCandidates->at(yIndex1);
		duplicate=false;
		yIndex2=0;
		while((!duplicate)&&(yIndex2<lenList->at(x))) {
			if(y==varList[x]->at(yIndex2)) {
				duplicate=true;
			}
			yIndex2++;
		}
		if(!duplicate) {
			varList[x]->push_back(y);
			lenList->at(x)++;
			total++;
		}
	}
}

void TVarListHandler::addToLine(int x, int yCandidate) {
	int yIndex;
	bool duplicate;
	duplicate=false;
	yIndex=0;
	while((!duplicate)&&(yIndex<lenList->at(x))) {
		if(yCandidate==varList[x]->at(yIndex)) {
			duplicate=true;
		}
		yIndex++;
	}
	if(!duplicate) {
		varList[x]->push_back(yCandidate);
		lenList->at(x)++;
		total++;
	}
}

void TVarListHandler::addToLine(int x, int yCandidate, bool testDuplicate) {
	if(testDuplicate) {
		addToLine(x,yCandidate);
		return;
	}
	
	varList[x]->push_back(yCandidate);
	lenList->at(x)++;
	total++;
}

void TVarListHandler::sort() {
	int x;
	for(x=0;x<res;x++) {
		sort(x);
	}
}

void TVarListHandler::sort(int x) {
	MSmergeSort<int>((int*) varList[x]->data(),(int) varList[x]->size(), &TVarListHandler::LowerEq);
}




bool TVarListHandler::LowerEq(int a, int b) {
	if(a<=b) {
		return true;
	}
	return false;
}


template <class T>
TVarListSignal<T>::TVarListSignal(TVarListHandler *_varList, T *_signal) {
	varList=_varList;
	signal=_signal;
	internalSignal=false;
}


template <class T>
TVarListSignal<T>::TVarListSignal(TVarListHandler *_varList, T init) {
	varList=_varList;
	signal=(T*) malloc(sizeof(T)*varList->total);
	for(int x=0;x<varList->total;x++) {
		signal[x]=init;
	}
	internalSignal=true;
}


template <class T>
TVarListSignal<T>::~TVarListSignal() {
	if (internalSignal) {
		if(signal!=NULL) {
			free(signal);
		}
	}
}


template <class T>
void TVarListSignal<T>::transcribeSorted(TVarListSignal<T> *src, T defaultValue) {
	int offset1,offset2,x,yIndex1,yIndex2,y1,y2;
	
	offset1=0;
	offset2=0;
	for(x=0;x<varList->res;x++) {
		yIndex1=0;
		yIndex2=0;
		while( ( yIndex1<varList->lenList->at(x) ) && ( yIndex2<src->varList->lenList->at(x) ) ) {
			y1=varList->varList[x]->at(yIndex1);
			y2=src->varList->varList[x]->at(yIndex2);
			if(y1==y2) {
				signal[offset1+yIndex1]=src->signal[offset2+yIndex2];
				yIndex1++;
				yIndex2++;
			} else {
				if(y1<y2) {
					signal[offset1+yIndex1]=defaultValue;
					yIndex1++;
				} else {
					yIndex2++;
				}
			}
		}
		while (yIndex1<varList->lenList->at(x)) {
			signal[offset1+yIndex1]=defaultValue;
			yIndex1++;			
		}
		offset1+=varList->lenList->at(x);
		offset2+=src->varList->lenList->at(x);
	}
	
}


template class TVarListSignal<int>;
template class TVarListSignal<double>;


