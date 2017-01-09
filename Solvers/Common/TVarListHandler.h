#ifndef TVARLISTHANDLER_H_
#define TVARLISTHANDLER_H_

#include<stdlib.h>
#include<vector>
#include<MergeSort.h>

using namespace std;

class TVarListHandler {
public:
	int res,total;
	vector<int> *lenList;
	vector<int> **varList;
	TVarListHandler();
	TVarListHandler(TVarListHandler *base);
	virtual ~TVarListHandler();
	void clear();
	void setupEmpty(int _res);
	void fillViaTranspose(TVarListHandler *transpose, int yres);
	void fillFromCSRIndexList(int *indices, int *indptr, int _res, int _total);
	void writeToCSRIndexList(int *indices, int *indptr);
	void addToLine(int x, vector<int>* yCandidates);
	void addToLine(int x, int yCandidate);
	void addToLine(int x, int yCandidate, bool testDuplicate);
	void sort();
	void sort(int x);
	
	static bool LowerEq(int a, int b);

};


template <class T>
class TVarListSignal {
public:
	TVarListHandler *varList;
	T *signal;
	bool internalSignal;
	
	TVarListSignal(TVarListHandler *_varList, T *_signal);
	TVarListSignal(TVarListHandler *_varList, T init);
	~TVarListSignal();
	
	void transcribeSorted(TVarListSignal<T> *src, T defaultValue);
};


#endif
