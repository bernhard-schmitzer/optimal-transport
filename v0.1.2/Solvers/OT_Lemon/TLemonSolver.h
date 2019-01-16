#ifndef TLemonSolver_H
#define TLemonSolver_H

#include<stdlib.h>
#include"TCouplingHandler.h"
#include<lemon/list_graph.h>
#include<lemon/network_simplex.h>
#include<lemon/cost_scaling.h>

using namespace lemon;
using namespace std;

typedef ListDigraph TGraph;
typedef TGraph::Node TNode;
typedef TGraph::Arc TArc;
typedef long int V;

typedef NetworkSimplex<TGraph,V,V> TAlgorithm;
typedef CostScaling<TGraph,V,V> TAlgorithmScaling;


template<class TCouplingHandlerType>
class TLemon_graph {
public:
	TGraph *G;
	TGraph::NodeMap<V> *supplyMap;
	TGraph::NodeMap<int> *idMap;
	vector<TNode> *xNodes,*yNodes;
	TGraph::ArcMap<int> *arcXMap, *arcYMap;
	TGraph::ArcMap<V> *cMap;
	int xres,yres;
	double *muX, *muY;
	double cScale;
	TCouplingHandlerType *CouplingHandler;

	TLemon_graph(TCouplingHandlerType *_CouplingHandler, double *_muX, double *_muY, double _cScale);
	~TLemon_graph();

	void setup();
	void setupNodes();
	void setupArcs();
	void setupArcs(bool reparametrizeCosts, double *alpha, double *beta);
	void deleteArcs();
};


class TLemonSolverBase {
public:

	static const int RSLT_ERROR_NOOPT=-101;
	static const int RSLT_ERROR_NOSOL=-102;
	static const int RSLT_ERROR_UNBOUNDED=-103;
	static const int RSLT_ERROR_UNDEFINED=-104;
	static const int NOT_IMPLEMENTED=-207;

	bool useDualVariables;

	double *muX, *muY;
	double *alpha, *beta;
	double objective;

	TLemonSolverBase(double *_muX, double *_muY);
	TLemonSolverBase(double *_muX, double *_muY, double *_alpha, double *_beta);

	virtual ~TLemonSolverBase() {};

	virtual int setup() { return NOT_IMPLEMENTED; };
	virtual int solve() { return NOT_IMPLEMENTED; };

	virtual int setupArcs() { return NOT_IMPLEMENTED; };
	virtual int setupArcs(__attribute__((unused)) bool dualOffset,
			__attribute__((unused)) double *alphaExt, __attribute__((unused)) double *betaExt)
			{ return NOT_IMPLEMENTED; };

	virtual int deleteArcs() { return NOT_IMPLEMENTED; };

	virtual int cleanup() { return NOT_IMPLEMENTED; };
	virtual double getObjective()  { return 0; };
	virtual double getObjective(__attribute__((unused)) bool dualOffset,
			__attribute__((unused)) double *alphaExt, __attribute__((unused)) double *betaExt)  { return 0; };
};

template <class TCouplingHandlerType, class TAlgorithmType>
class TLemonSolver : public TLemonSolverBase {
public:

	TCouplingHandlerType *CouplingHandler;
	TLemon_graph<TCouplingHandlerType> *graph;
	double cScale;


	TLemonSolver(TCouplingHandlerType *_CouplingHandler, double *_muX, double *_muY, double _cScale);
	TLemonSolver(TCouplingHandlerType *_CouplingHandler, double *_muX, double *_muY, double *_alpha, double *_beta, double _cScale);
	~TLemonSolver();

	int setup();
	int setupArcs();
	int setupArcs(bool dualOffset, double *alphaExt, double *betaExt);
	int deleteArcs();

	int solve();
	double getObjective();
	double getObjective(bool dualOffset, double *alphaExt, double *betaExt);

};


///////////////////////////////////////////////////////////////////////////////////
// Methods TLemon_graph<TCouplingHandlerType>

template <class TCouplingHandlerType>
TLemon_graph<TCouplingHandlerType>::TLemon_graph(TCouplingHandlerType *_CouplingHandler, double *_muX, double *_muY,
		double _cScale) {
	CouplingHandler=_CouplingHandler;
	xres=CouplingHandler->xres;
	yres=CouplingHandler->yres;
	muX=_muX;
	muY=_muY;
	cScale=_cScale;

	G=NULL;
	supplyMap=NULL;
	idMap=NULL;
	xNodes=NULL;
	yNodes=NULL;
	arcXMap=NULL;
	arcYMap=NULL;
	cMap=NULL;
}


template <class TCouplingHandlerType>
TLemon_graph<TCouplingHandlerType>::~TLemon_graph() {
	if (cMap!=NULL) {
		delete cMap;
	}
	if (arcXMap!=NULL) {
		delete arcXMap;
	}
	if (arcYMap!=NULL) {
		delete arcYMap;
	}
	if (xNodes!=NULL) {
		delete xNodes;
	}
	if (yNodes!=NULL) {
		delete yNodes;
	}
	if (idMap!=NULL) {
		delete idMap;
	}
	if (supplyMap!=NULL) {
		delete supplyMap;
	}
	if (G!=NULL) {
		delete G;
	}
}

template <class TCouplingHandlerType>
void TLemon_graph<TCouplingHandlerType>::setup() {
		setupNodes();
		setupArcs();
}

template <class TCouplingHandlerType>
void TLemon_graph<TCouplingHandlerType>::setupNodes() {
	int xpos,ypos;

	G=new TGraph();
	G->reserveNode(xres+yres);
	supplyMap= new TGraph::NodeMap<V>(*G);
	idMap= new TGraph::NodeMap<int>(*G);
	xNodes = new vector<TNode>(xres);
	yNodes = new vector<TNode>(yres);

	arcXMap=new TGraph::ArcMap<int>(*G);
	arcYMap=new TGraph::ArcMap<int>(*G);
	cMap=new TGraph::ArcMap<V>(*G);

	for(xpos=0;xpos<xres;xpos++) {
		(*xNodes)[xpos]=G->addNode();
		(*supplyMap)[(*xNodes)[xpos]]=(V) muX[xpos];
		(*idMap)[(*xNodes)[xpos]]=xpos;
	}
	for(ypos=0;ypos<yres;ypos++) {
		(*yNodes)[ypos]=G->addNode();
		(*supplyMap)[(*yNodes)[ypos]]=-((V) muY[ypos]);
		(*idMap)[(*yNodes)[ypos]]=ypos;
	}
}

template <class TCouplingHandlerType>
void TLemon_graph<TCouplingHandlerType>::setupArcs() {
	setupArcs(false,NULL,NULL);
}

template <class TCouplingHandlerType>
void TLemon_graph<TCouplingHandlerType>::setupArcs(bool reparametrizeCosts, double *alpha, double *beta) {
	int xpos,rowlen,yIndex,ypos;
	TArc arc;

	G->reserveArc(CouplingHandler->total);


	// iterate through all rows
	for(xpos=0;xpos<xres;xpos++) {
		// for each row, iterate through columns of that row
		rowlen=CouplingHandler->getRowLength(xpos);
		// set xside, yside and cost function for all arcs
		for(yIndex=0;yIndex<rowlen;yIndex++) {
			ypos=CouplingHandler->getColNr(xpos,yIndex);
			arc=G->addArc((*xNodes)[xpos],(*yNodes)[ypos]);
			if (reparametrizeCosts) {
				(*cMap)[arc]=(V)((CouplingHandler->getCRow(xpos,yIndex)-alpha[xpos]-beta[ypos])/cScale);
			} else {
				(*cMap)[arc]=(V)(CouplingHandler->getCRow(xpos,yIndex)/cScale);
			}
			(*arcXMap)[arc]=xpos;
			(*arcYMap)[arc]=yIndex;
		}
	}
}


template <class TCouplingHandlerType>
void TLemon_graph<TCouplingHandlerType>::deleteArcs() {
	for(TGraph::ArcIt a(*G);a!=INVALID;++a) {
		G->erase(a);
	}
}


////////////////////////////////////////////////////////////////////////////
// TLemonSolver<TCouplingHandlerType> Methods

template<class TCouplingHandlerType, class TAlgorithmType>
TLemonSolver<TCouplingHandlerType, TAlgorithmType>::TLemonSolver(
		TCouplingHandlerType* _CouplingHandler, double* _muX, double* _muY, double _cScale) : TLemonSolverBase(_muX,_muY) {
	CouplingHandler=_CouplingHandler;
	graph=NULL;
	cScale=_cScale;

}

template<class TCouplingHandlerType, class TAlgorithmType>
TLemonSolver<TCouplingHandlerType, TAlgorithmType>::TLemonSolver(
		TCouplingHandlerType* _CouplingHandler, double* _muX, double* _muY,
		double* _alpha, double* _beta, double _cScale) : TLemonSolverBase(_muX,_muY,_alpha,_beta) {
	CouplingHandler=_CouplingHandler;
	graph=NULL;
	cScale=_cScale;
}

template<class TCouplingHandlerType, class TAlgorithmType>
int TLemonSolver<TCouplingHandlerType, TAlgorithmType>::setup() {
	graph=new TLemon_graph<TCouplingHandlerType>(CouplingHandler,muX,muY,cScale);
	graph->setup();
	return 0;
}


template<class TCouplingHandlerType, class TAlgorithmType>
int TLemonSolver<TCouplingHandlerType, TAlgorithmType>::setupArcs() {
	graph->setupArcs();
	return 0;
}

template<class TCouplingHandlerType, class TAlgorithmType>
int TLemonSolver<TCouplingHandlerType, TAlgorithmType>::setupArcs(bool dualOffset, double *alphaExt, double *betaExt) {
	graph->setupArcs(dualOffset, alphaExt, betaExt);
	return 0;
}


template<class TCouplingHandlerType, class TAlgorithmType>
int TLemonSolver<TCouplingHandlerType, TAlgorithmType>::deleteArcs() {
	graph->deleteArcs();
	return 0;
}


template<class TCouplingHandlerType, class TAlgorithmType>
int TLemonSolver<TCouplingHandlerType, TAlgorithmType>::solve() {
		TAlgorithmType algorithm(*graph->G);
		algorithm.supplyMap(*graph->supplyMap);
		algorithm.costMap(*graph->cMap);

		typename TAlgorithmType::ProblemType status;
		status=algorithm.run();

		switch(status) {
		case TAlgorithm::OPTIMAL:
			break;
		case TAlgorithm::INFEASIBLE:
			return RSLT_ERROR_NOSOL;
			break;
		case TAlgorithm::UNBOUNDED:
			return RSLT_ERROR_UNBOUNDED;
			break;
		default:
			return RSLT_ERROR_UNDEFINED;
		}

		int xpos,yIndex;
		for(TGraph::ArcIt a(*graph->G);a!=INVALID;++a) {
			xpos=(*graph->arcXMap)[a];
			yIndex=(*graph->arcYMap)[a];
			CouplingHandler->setMuRow(xpos,yIndex,(double) algorithm.flow(a));
		}

		if (useDualVariables) {
			int i;
			for(i=0;i<graph->xres;i++) {
				alpha[i]=-((double) algorithm.potential(graph->xNodes->at(i)))*cScale;
			}
			for(i=0;i<graph->yres;i++) {
				beta[i]=((double) algorithm.potential(graph->yNodes->at(i)))*cScale;
			}
		}

		objective=algorithm.totalCost()*cScale;


		return 0;
	}

template<class TCouplingHandlerType, class TAlgorithmType>
TLemonSolver<TCouplingHandlerType, TAlgorithmType>::~TLemonSolver() {
	if (graph!=NULL) {
		delete graph;
	}
}

template<class TCouplingHandlerType, class TAlgorithmType>
double TLemonSolver<TCouplingHandlerType, TAlgorithmType>::getObjective() {
	return objective;
}

template<class TCouplingHandlerType, class TAlgorithmType>
double TLemonSolver<TCouplingHandlerType, TAlgorithmType>::getObjective(bool dualOffset, double *alphaExt, double *betaExt) {
	if (dualOffset) {
		double realObjective;
		int i;
		realObjective=objective;
		for(i=0;i<CouplingHandler->xres;i++) {
			realObjective+=alphaExt[i]*muX[i];
		}
		for(i=0;i<CouplingHandler->yres;i++) {
			realObjective+=betaExt[i]*muY[i];
		}
		return realObjective;
	}
	return objective;
}

#endif
