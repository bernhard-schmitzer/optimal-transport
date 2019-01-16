#ifndef TSINKHORN_H_
#define TSINKHORN_H_

#include<stdlib.h>
#include"TMultiVarListHandler.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// different data structures for potentially handling multi-marginal kernels
// the further classes (hierarchical constraint search, cost function evaluation, cell refinement, etc.)
// can then be instantiated with different data structures
// the interface is given by
// addVariables, getIterator, iterate_pos, iterate_both

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Multi_CSR: one var list for each marginal, in CSR style, but now each col-position is no longer just an integer
// but an array of integers, one for all other marginals
template<class TMultiVarListHandlerType, class T>
class TKernelHandler_MultiCSR_Base {
public:
	struct iterator {
		int row, col;
	};
	TMultiVarListHandlerType **varLists;
	int dim;

	TKernelHandler_MultiCSR_Base(int _dim, int *_res);
	~TKernelHandler_MultiCSR_Base();

	void addVariables(int *pos, T value);

	inline iterator getIterator() {
		iterator result;
		result.col=0;
		result.row=0;
		return result;
	}
	inline void step_iterator(iterator *it) {
		if(it->row<varLists[0]->res) {
			if(it->col < varLists[0]->lenList->at(it->row)-1) {
				it->col++;
			} else {
				it->row++;
				it->col=0;
			}
		}
	}
	inline bool check_iterator(iterator *it) {
		if(it->row<varLists[0]->res) {
			if(it->col < varLists[0]->lenList->at(it->row)) {
				return true;
			}
		}
		return false;
	}
	inline bool iterate_pos(iterator *it,int * const pos) {
		if(! check_iterator(it)) { return false; };
		// copy current pos into pos
		pos[0]=it->row;
		copy(varLists[0]->varList[it->row]->at(it->col),varLists[0]->varList[it->row]->at(it->col)+dim-1,pos+1);
		// increase iterator
		step_iterator(it);
		return true;
	}

	inline bool iterate_both(iterator *it, int * const pos, T * const sig) {
		if(! check_iterator(it)) { return false; };
		// copy current pos into pos
		pos[0]=it->row;
		copy(varLists[0]->varList[it->row]->at(it->col),varLists[0]->varList[it->row]->at(it->col)+dim-1,pos+1);
		// assign signal
		sig[0]=varLists[0]->signalList[it->row]->at(it->col);
		step_iterator(it);
		return true;
	}

};

//typedef TKernelHandler_MultiCSR_Base<TMultiVarListHandler<double>, double >
//	TKernelHandler_MultiCSR;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// standard 2-marginal CSR kernel handler.
// basically just a wrapper around the old TVarListHandler for the interface methods
// does not come with a signal handler. is just used to store indices (e.g. for refinement, and re-evaluation)
class TKernelHandler_CSR {
public:
	struct iterator {
		int row, col;
	};
	TVarListHandler *varList;
	int dim;

	TKernelHandler_CSR(int _dim, int *_res);
	TKernelHandler_CSR(TVarListHandler *_varList);
	~TKernelHandler_CSR();

	void addVariables(int *pos, __attribute__((unused))  double value);

	inline iterator getIterator() {
		iterator result;
		result.col=0;
		result.row=0;
		return result;
	}
	inline void step_iterator(iterator *it) {
		if(it->row<varList->res) {
			if(it->col < varList->lenList->at(it->row)-1) {
				it->col++;
			} else {
				it->row++;
				it->col=0;
			}
		}
	}
	inline bool check_iterator(iterator *it) {
		if(it->row<varList->res) {
			if(it->col < varList->lenList->at(it->row)) {
				return true;
			}
		}
		return false;
	}
	inline bool iterate_pos(iterator *it,int * const pos) {
		if(! check_iterator(it)) { return false; };
		// copy current pos into pos
		pos[0]=it->row;
		pos[1]=varList->varList[it->row]->at(it->col);
		// increase iterator
		step_iterator(it);
		return true;
	}


};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// pos based multi-marginal data structure. for each kernel entry a value and its full dimension (list of ints)
// is stored

template<class T>
class TKernelHandler_Pos_Base {
public:
	vector<int*> *varList;
	vector<T> *signalList;
	int dim;
	int total;

	struct iterator {
		int pos;
	};

	TKernelHandler_Pos_Base(int _dim, __attribute__((unused)) int *_res) {
		dim=_dim;
		total=0;
		varList=new vector<int*>();
		signalList=new vector<T>();

	}
	~TKernelHandler_Pos_Base() {
		int i;
		for(i=0;i<total;i++) {
			free(varList->at(i));
		}
		delete(varList);
		delete(signalList);
	}


	inline void addVariables(int *pos, T value) {
		int *poscopy=(int*) malloc(sizeof(int)*dim);
		copy(pos,pos+dim,poscopy);
		varList->push_back(poscopy);
		signalList->push_back(value);
		total++;
	}
	inline iterator getIterator() {
		iterator result;
		result.pos=0;
		return result;
	}
	inline void step_iterator(iterator *it) {
		it->pos++;
	}
	inline bool check_iterator(iterator *it) {
		if(it->pos<((int) varList->size())) {
			return true;
		}
		return false;
	}
	inline bool iterate_pos(iterator *it,int * const pos) {
		if(! check_iterator(it)) { return false; };
		copy(varList->at(it->pos),varList->at(it->pos)+dim,pos);
		step_iterator(it);
		return true;
	}
	inline bool iterate_both(iterator *it,int * const pos, T * const sig) {
		if(! check_iterator(it)) { return false; };
		copy(varList->at(it->pos),varList->at(it->pos)+dim,pos);
		sig[0]=signalList->at(it->pos);
		step_iterator(it);
		return true;
	}

	void writeToPosIndexList(T *data, int *pos) {
		int x,d;
		for(x=0;x<total;x++) {
			data[x]=signalList->at(x);
			for(d=0;d<dim;d++) {
				pos[dim*x+d]=varList->at(x)[d];
			}
		}
	}

	void fillFromPosIndexList(T *data, int *pos, int _total) {
		int x,d;

		// possibly free old var list stuff
		for(x=0;x<total;x++) {
			free(varList->at(x));
		}

		total=_total;
		varList->resize(total);
		signalList->resize(total);

		for(x=0;x<total;x++) {
			signalList->at(x)=data[x];
			varList->at(x)=(int*) malloc(sizeof(int)*dim);
			for(d=0;d<dim;d++) {
				varList->at(x)[d]=pos[dim*x+d];
			}
		}
	}

	void rescale(T **scalings) {
		int x,d;
		for(x=0;x<total;x++) {
			for(d=0;d<dim;d++) {
				signalList->at(x)*=scalings[d][varList->at(x)[d]];
			}
		}
	}

	void exponentiate(T eps) {
		int x;
		for(x=0;x<total;x++) {
			signalList->at(x)=exp(-signalList->at(x)/eps);
		}
	}

};

typedef TKernelHandler_Pos_Base<double >
	TKernelHandler_Pos;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// THierarchicalDualConstraintSearch
/* detects relevant (multi-marginal)-kernel entries, based on an effective cost function (cost + potentials),
    and constructs the corresponding varList
	relevance is tested with a scalar threshold (slack)
	
	the search is implemented as follows:
		initialization is done in SearchDualConstraints(...), which also launches the top-level search by calling
			SearchDualConstraints_checkCell(...)
		SearchDualConstraints_checkCell(...) checks the slack of a given cell on a given level. if the cell is relevant:
			if on the finest level, add the relevant cell to the varLists
			if not on finest level, call SearchDualConstraints_checkChildren(...), which in turn calls
				SearcHDualConstraints_checkCell(...) on all children of this cell
	*/


// derived class, specialized on CostFunctionProvider.
template<class TKernelHandlerType, class T>
class THierarchicalDualConstraintSearchBase {
public:
	int dim; // number of marginals
	THierarchicalPartition **partitions; // partitions over marginals

	TMultiCostFunctionProvider *costProvider;
	
	THierarchicalDualConstraintSearchBase(int _dim, THierarchicalPartition **_partitions,
			TMultiCostFunctionProvider *_costProvider);
	
	TKernelHandlerType* SearchDualConstraints(T const slack, int layerBottom);
	void SearchDualConstraints_checkCell(
			T const slack, TKernelHandlerType * const kernel, int layerBottom,
			int layer, int *pos);
	void SearchDualConstraints_checkChildren(
			T const slack, TKernelHandlerType * const kernel,
			int layerBottom, int layerParent, int *posParent);
};

typedef THierarchicalDualConstraintSearchBase<TKernelHandler_Pos, double >
		THierarchicalDualConstraintSearch_Pos;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TMultiVarListRefinement
/* see .cpp file at Refine(...) for a description . */

template<class TKernelHandlerTypeIn, class TKernelHandlerTypeOut, class T>
class TMultiVarListRefinementBase {
public:

	int dim;
	THierarchicalPartition **partitions;
	TMultiCostFunctionProvider *costProvider;

	TMultiVarListRefinementBase(int _dim, THierarchicalPartition **_partitions,
			TMultiCostFunctionProvider *_costProvider);

	TKernelHandlerTypeOut* Refine(TKernelHandlerTypeIn *kernelTop, int layerTop);
	void Refine_refineCell(TKernelHandlerTypeOut * const kernel, int layerTop, int *posParent);
};

typedef TMultiVarListRefinementBase<TKernelHandler_CSR,TKernelHandler_Pos, double >
		TMultiVarListRefinement_CSR_Pos;

typedef TMultiVarListRefinementBase<TKernelHandler_Pos,TKernelHandler_Pos, double >
		TMultiVarListRefinement_Pos_Pos;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TMultiVarListEvaluation

template<class TKernelHandlerTypeIn, class TKernelHandlerTypeOut, class T>
class TMultiVarListEvaluationBase {
public:

	int dim;
	THierarchicalPartition **partitions;
	TMultiCostFunctionProvider *costProvider;

	TMultiVarListEvaluationBase(int _dim, THierarchicalPartition **_partitions,
			TMultiCostFunctionProvider *_costProvider);

	TKernelHandlerTypeOut* Evaluate(TKernelHandlerTypeIn *kernelOld, int layerTop);
};

typedef TMultiVarListEvaluationBase<TKernelHandler_CSR, TKernelHandler_Pos, double >
	TMultiVarListEvaluation_CSR_Pos;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TKLProjection
/* implements a convolution of a hyper-kernel with a number of dual potentials and the corresponding
	KL-projection onto the marginal constraint.
	
	scaling contains list of scalings. scaling[0] is the scaling of the space on which k is based. scaling[1] -> scaling[dim-1] 
	contains the other potentials */


template<class TKernelHandlerType, class T>
class TKLProjectionBase {
public:
	int dim;
	TKernelHandlerType *k;
	T **scaling;
	int *res;

	TKLProjectionBase(TKernelHandlerType *_k, T **_scaling, int *_res);

	void getConvolution(int axis, double *muOut);
	void project(int axis, double *mu);
	void projectSoft_KL(int axis, double *mu, double q);
	void getMarginal(int axis, double *muOut);
};

typedef TKLProjectionBase<TKernelHandler_Pos, double > TKLProjection_Pos;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TDenseKernel
/* constructs a set of dense kernels in the TMultiVarList representation. useful for starting the Sinkhorn
 * algorithm on the coarsest scale.
 */

template<class TKernelHandlerType, class T>
class TDenseCostGetterBase {
public:
	int dim;
	int layer;
	THierarchicalPartition **partitions;
	TMultiCostFunctionProvider *costProvider;

	TDenseCostGetterBase(int _dim, THierarchicalPartition **_partitions, TMultiCostFunctionProvider *_costProvider,
			int _layer);

	TKernelHandlerType* getCost();
};

typedef TDenseCostGetterBase<TKernelHandler_Pos, double > TDenseCostGetter_Pos;
//typedef TDenseCostGetterBase<TKernelHandler_CSR, double > TDenseCostGetter_CSR;


#endif
