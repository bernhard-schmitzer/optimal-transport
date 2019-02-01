#include "TVarListHandler.h"

TVarListHandler::TVarListHandler() {
	res=0;
	total=0;
}

// copy constructor
TVarListHandler::TVarListHandler(const TVarListHandler* const base) : lenList(base->lenList), varList(base->varList) {
	res=base->res;
	total=base->total;
}


// initialize with CSR Index List
TVarListHandler::TVarListHandler(int *indices, int *indptr, int _res, int _total, bool sorted) {

	fillFromCSRIndexList(indices,indptr,_res,_total);
	if (sorted) {
		sort();
	}
	
}

TVarListHandler::~TVarListHandler() {
}

void TVarListHandler::clear() {
	lenList.resize(0);
	varList.resize(0);
	res=0;
	total=0;
}

void TVarListHandler::setupEmpty(int _res) {
	res=_res;
	lenList.resize(res);
	varList.resize(res);
}

void TVarListHandler::fillViaTranspose(const TVarListHandler* const transpose, int yres) {
	setupEmpty(yres);

	for(int x=0;x<transpose->res;x++) {
		for(int yIndex=0;yIndex<transpose->lenList[x];yIndex++) {
			int y=transpose->varList[x][yIndex];
			varList[y].push_back(x);
		}
	}

	total=0;
	for(int y=0;y<yres;y++) {
		lenList[y]=(int) varList[y].size();
		total+=lenList[y];
	}

}

void TVarListHandler::fillFromCSRIndexList(int *indices, int *indptr, int _res, int _total) {
	setupEmpty(_res);
	total=_total;
	for(int x=0;x<_res;x++) {
		int rowLen,offset;
		
		rowLen=indptr[x+1]-indptr[x];
		offset=indptr[x];
		lenList[x]=rowLen;
		varList[x].resize(rowLen);
		for(int y=0;y<rowLen;y++) {
			varList[x][y]=indices[offset+y];
		}
	}
}

void TVarListHandler::writeToCSRIndexList(int *indices, int *indptr) {
	// write varList data to CSR (as in scipy.sparse.csr_matrix) list. indices must have size total, indptr must have size res+1
	int x,offset;
	indptr[0]=0;
	offset=0;
	for(x=0;x<res;x++) {
		for(int yIndex=0;yIndex<lenList[x];yIndex++) {
			indices[offset]=varList[x][yIndex];
			offset++;
		}
		indptr[x+1]=offset;
	}
	
}

int TVarListHandler::merge(TVarListHandler *addition) {
	if(res!=addition->res) {
		// res of two var lists must coincide. otherwise indicate error
		return ERR_BASE_VARLIST_MERGE_RESFAIL;
	}
	TVarListHandler::TIterator it=addition->iterationInitialize();
	while(addition->iterate(&it)) {
		addToLine(it.x,it.y,true);
	}	
	
	return 0;
}

int TVarListHandler::mergeSelected(TVarListHandler *addition, bool *selection) {
	if(res!=addition->res) {
		// res of two var lists must coincide. otherwise indicate error
		return ERR_BASE_VARLIST_MERGE_RESFAIL;
	}
	TVarListHandler::TIterator it=addition->iterationInitialize();
	while(addition->iterate(&it)) {
		if(selection[it.offset]) {
			addToLine(it.x,it.y,true);
		}
	}	

	return 0;
	
}


void TVarListHandler::addToLine(int x, const std::vector<int>& yCandidates) {
	for(std::size_t yIndex1=0;yIndex1< yCandidates.size();yIndex1++) {
		int y=yCandidates[yIndex1];
		bool duplicate=false;
		int yIndex2=0;
		while((!duplicate)&&(yIndex2<lenList[x])) {
			if(y==varList[x][yIndex2]) {
				duplicate=true;
			}
			yIndex2++;
		}
		if(!duplicate) {
			varList[x].push_back(y);
			lenList[x]++;
			total++;
		}
	}
}

void TVarListHandler::addToLine(int x, const int* const yCandidates, int nCandidates) {
	for(int yIndex1=0;yIndex1< nCandidates;yIndex1++) {
		int y=yCandidates[yIndex1];
		bool duplicate=false;
		int yIndex2=0;
		while((!duplicate)&&(yIndex2<lenList[x])) {
			if(y==varList[x][yIndex2]) {
				duplicate=true;
			}
			yIndex2++;
		}
		if(!duplicate) {
			varList[x].push_back(y);
			lenList[x]++;
			total++;
		}
	}
}


void TVarListHandler::addToLine(int x, int yCandidate) {
	bool duplicate=false;
	int yIndex=0;
	while((!duplicate)&&(yIndex<lenList[x])) {
		if(yCandidate==varList[x][yIndex]) {
			duplicate=true;
		}
		yIndex++;
	}
	if(!duplicate) {
		varList[x].push_back(yCandidate);
		lenList[x]++;
		total++;
	}
}

void TVarListHandler::addToLine(int x, int yCandidate, bool testDuplicate) {
	if(testDuplicate) {
		addToLine(x,yCandidate);
		return;
	}
	
	varList[x].push_back(yCandidate);
	lenList[x]++;
	total++;
}

void TVarListHandler::sort() {
	int x;
	for(x=0;x<res;x++) {
		sort(x);
	}
}

void TVarListHandler::sort(int x) {
	//MSmergeSort<int>((int*) varList[x]->data(),(int) varList[x]->size(), &TVarListHandler::LowerEq);
	std::sort(varList[x].data(),varList[x].data()+varList[x].size(), &TVarListHandler::LowerEq);
}




bool TVarListHandler::LowerEq(int a, int b) {
	if(a<=b) {
		return true;
	}
	return false;
}


// iteration

TVarListHandler::TIterator TVarListHandler::iterationInitialize() {
	return {0,0,0,0,true};
}
bool TVarListHandler::iterate(TVarListHandler::TIterator *it) {
	// if no rows, then always terminate iteration immediately
	if(res==0) {
		return false;
	}

	if(!it->iterationInitialized) {
		// if not first call after initialization, increase col index
		it->yIndex++;
		// also increase global variable nr index
		it->offset++;
	} else {
		it->iterationInitialized=false;
	}
	if(it->yIndex>=lenList[it->x]) {
		// if end of row is reached, increase row index
		it->x++;
		if(it->x>=res) {
			// if number of rows is exceeded, return false
			return false;
		}
		while(lenList[it->x]==0) {
			// until non-empty row is found, further increase row index
			it->x++;
			if(it->x>=res) {
				return false;
			}
		}
		// if non-empty row is found, reset col index
		it->yIndex=0;
	}
	// extract absolute col index
	it->y=varList[it->x][it->yIndex];
	return true;	
	
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <class T>
TVarListSignal<T>::TVarListSignal(TVarListHandler *_varList, T *_signal) {
	varList=_varList;
	signal=_signal;
	internalSignal=false;
	
	offsets=NULL;
	computedOffsets=false;
}


template <class T>
TVarListSignal<T>::TVarListSignal(TVarListHandler *_varList, T init) {
	varList=_varList;
	signal=(T*) malloc(sizeof(T)*varList->total);
	for(int x=0;x<varList->total;x++) {
		signal[x]=init;
	}
	internalSignal=true;
	offsets=NULL;
	computedOffsets=false;
}


template <class T>
TVarListSignal<T>::~TVarListSignal() {
	if (internalSignal) {
		if(signal!=NULL) {
			free(signal);
		}
	}
	if (computedOffsets) {
		if(offsets!=NULL) {
			free(offsets);
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
		while( ( yIndex1<varList->lenList[x] ) && ( yIndex2<src->varList->lenList[x] ) ) {
			y1=varList->varList[x][yIndex1];
			y2=src->varList->varList[x][yIndex2];
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
		while (yIndex1<varList->lenList[x]) {
			signal[offset1+yIndex1]=defaultValue;
			yIndex1++;			
		}
		offset1+=varList->lenList[x];
		offset2+=src->varList->lenList[x];
	}
	
}

template <class T>
void TVarListSignal<T>::computeOffsets() {
	
	computedOffsets=true;

	offsets=(int*) malloc(sizeof(int)*varList->total);
	offsets[0]=0;
	for(int i=0;i<varList->res-1;i++) {
		offsets[i+1]=offsets[i]+varList->lenList[i];
	}
}

template <class T>
void TVarListSignal<T>::write(int x, int y, T value) {
	int yIndex=0;
	while(true) {
		if (yIndex>=varList->lenList[x]) {
			eprintf("ERROR: TVarListSignal::write failed because y element was not found.\n");
			return;
		}
		if (varList->varList[x][yIndex]==y) {
			signal[offsets[x]+yIndex]=value;
			return;
		}
		yIndex++;
	}
}

template <class T>
void TVarListSignal<T>::writeFromTranspose(TVarListSignal<T> *transpose) {
	TVarListHandler::TIterator it=transpose->varList->iterationInitialize();
	while(transpose->varList->iterate(&it)) {
		write(it.y,it.x,transpose->signal[it.offset]);
	}	
	
}


template class TVarListSignal<bool>;
template class TVarListSignal<int>;
template class TVarListSignal<double>;


TVarListHandler* GetFullVarList(int xres, int yres) {
	TVarListHandler *result;
	result=new TVarListHandler();
	result->setupEmpty(xres);
	for(int i=0;i<xres;i++) {
		result->lenList[i]=yres;
		result->varList[i].resize(yres);
		for(int j=0;j<yres;j++) {
			result->varList[i][j]=j;
		}
	}
	result->total=xres*yres;
	return result;
}

bool VarListTools_HasEmptyRows(TVarListHandler *vars) {
	for(int i=0;i<vars->res;i++) {
		if(vars->lenList[i]==0) {
			return true;
		}
	}
	return false;
}
