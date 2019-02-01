#ifndef __Tools_H
#define __Tools_H

#include<cstdlib>
#include<fstream>
#include<iterator>
#include<vector>
#include<cmath>
#include<cstring>
#include<forward_list>
#include<algorithm>

#include<Common/ErrorCodes.h>
#include<Common/PythonTypes.h>
//#include<Common/Verbose.h>

void doubleArrayCopy(double *a, double *b, int n);
void doubleArrayScale(double *a, double b, int n);

double doubleArrayMin(double *data, int res);
double doubleArrayMax(double *data, int res);

int intArrayMin(int *data, int res);
int intArrayMax(int *data, int res);

void freeTDoubleMatrix(TDoubleMatrix *a);

// some basic data structure structs

struct TSparseCSRContainer {
	std::vector<double> data;
	std::vector<int> indices;
	std::vector<int> indptr;
	int xres,yres,nonZeros;
};

struct TSparsePosContainer {
	std::vector<int> posStart;
	std::vector<int> posEnd;
	std::vector<double> mass;
	int xres,yres,nParticles;
};

struct TParticleContainer {
	std::vector<double> pos;
	std::vector<double> mass;
	int nParticles;
	int dim;
};


////////////////////////////////

// a few geometric aux functions that are often used
double EUCL_innerProduct(const double * const a, const double * const b, const int n);
void EUCL_lincomb(const double * const a, const double * const b, double * const c, const double sa, const double sb, const int n);
double EUCL_lincombSqr(const double * const a, const double * const b, const double sa, const double sb, const int n);
double EUCL_len(const double * const a, const int n);
///////////////////////////////////////////////////////
// read raw double data from file
// (used in examples)

//int readDoubleFile(const char* filename, double **data, int *len);
template <typename V>
std::vector<V> readFile(const char* filename);
template <typename V>
int writeFile(const char* filename, const V * const data, const int len);

////////////////////////////////////////////////////////////////////////////////////////////
// list template for hierarchical search algorithms

template <class TElementType>
class THierarchicalSearchList {
public:
	int size;
	std::forward_list<TElementType> data;
	
	THierarchicalSearchList(TElementType *newData, int n);
	~THierarchicalSearchList();
	void merge(TElementType *newData, int n);
};

template <class TElementType>
THierarchicalSearchList<TElementType>::THierarchicalSearchList(TElementType *newData, int n) {
	// construct linked list from array data
	
	// first sort
	std::sort(newData,newData+n);
	
	// then iteratively add elements
	typename std::forward_list<TElementType>::iterator it=data.before_begin();
	for(int i=0;i<n;i++) {
		it=data.insert_after(it,newData[i]);
	}
	// set number of elements
	size=n;
}

template <class TElementType>
THierarchicalSearchList<TElementType>::~THierarchicalSearchList() {}

template <class TElementType>
void THierarchicalSearchList<TElementType>::merge(TElementType *newData, int n) {
	// first sort
	std::sort(newData,newData+n);

	// if current list is empty, simply add new elements
	if(size==0) {
		// then iteratively add elements
		typename std::forward_list<TElementType>::iterator it=data.before_begin();
		for(int i=0;i<n;i++) {
			it=data.insert_after(it,newData[i]);
		}
		// set size
		size=n;
		return;
	}
	
	// keep two iterators pointing to subsequent elements
	typename std::forward_list<TElementType>::iterator it=data.before_begin();
	typename std::forward_list<TElementType>::iterator itNext=data.begin();

	int i=0;
	// go through list, as long as new elements remain to be added
	while((itNext!=data.end()) && (i<n)) {
		if(newData[i]<*itNext) {
			// if new element should be between it and itNext
			
			// insert element after it, set it to new element
			it=data.insert_after(it,newData[i]);
			// set itNext to subsequent element
			itNext=it;
			itNext++;
			// focus on next new element
			i++;
		} else {
			// otherwise, simply go through list
			it++;
			itNext++;
		}
	}
	while(i<n) {
		// if some new elements are still left, add them to end of list
		it=data.insert_after(it,newData[i]);
		i++;
	}
	// set new size
	size+=n;
}

#endif
