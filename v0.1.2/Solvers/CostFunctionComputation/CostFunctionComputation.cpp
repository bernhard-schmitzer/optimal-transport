#include<stdlib.h>
#include<PythonTypes.h>
#include"TVarListHandler.h"

#include<THierarchicalPartition.h>
//#include<THierarchicalPartitionInterface.h>

using namespace std;

const int MODE_MIN=0;
const int MODE_MAX=1;

extern "C" {

int refineVarListRequestMemory(long int PartitionXAddr, long int PartitionYAddr,
		TInteger32Matrix *VLXindices, TInteger32Matrix *VLXindptr, int startLayerX, int startLayerY,
		TInteger32Matrix *output) {
	// compute how many variables will be in the refined var list

	THierarchicalPartition *partitionX=(THierarchicalPartition*) PartitionXAddr;
	THierarchicalPartition *partitionY=(THierarchicalPartition*) PartitionYAddr;

	int result,resultRow;
	int x, yIndex, y;
	int xres;

	xres=partitionX->layers[startLayerX]->nCells;

	result=0;
	// go through all rows in the partition
	for(x=0;x<xres;x++) {
		// in each row go through all cols
		resultRow=0;
		// go through yIndices in CSR col
		for(yIndex=VLXindptr->data[x];yIndex<VLXindptr->data[x+1];yIndex++) {
			// look up y position
			y=VLXindices->data[yIndex];
			// add number of children in that y cell
			resultRow+=partitionY->layers[startLayerY]->nChildren[y];
		}
		// multiply this result by number of children in x cell, add to total
		result+=resultRow*partitionX->layers[startLayerX]->nChildren[x];
	}

	output->data[0]=partitionX->layers[startLayerX+1]->nCells;
	output->data[1]=result;
	return 0;
}


int refineVarList(long int PartitionXAddr, long int PartitionYAddr,
		TInteger32Matrix *VLXindices, TInteger32Matrix *VLXindptr, int startLayerX, int startLayerY,
		TInteger32Matrix *VLXindicesOut, TInteger32Matrix *VLXindptrOut) {
	// compute how many variables will be in the refined var list

	THierarchicalPartition *partitionX=(THierarchicalPartition*) PartitionXAddr;
	THierarchicalPartition *partitionY=(THierarchicalPartition*) PartitionYAddr;

	TPartitionLayer *layerY, *layerXFine;
	layerY=partitionY->layers[startLayerY];
	layerXFine=partitionX->layers[startLayerX+1];

	int offsetInd,offsetSep;
	int x, yIndex, y,yFineIndex,xFine,yFine;
	int xresFine;
	xresFine=layerXFine->nCells;

	VLXindptrOut->data[0]=0;
	offsetInd=0;
	offsetSep=1;
	// go through all rows in the partition
	for(xFine=0;xFine<xresFine;xFine++) {
		// for each xFine determine parent
		x=layerXFine->parent[xFine];
		// in parent row go through all cols
		// go through yIndices in CSR col
		for(yIndex=VLXindptr->data[x];yIndex<VLXindptr->data[x+1];yIndex++) {
			// look up y position
			y=VLXindices->data[yIndex];
			// go through children of y
			for(yFineIndex=0;yFineIndex<layerY->nChildren[y];yFineIndex++) {
				yFine=layerY->children[y][yFineIndex];
				VLXindicesOut->data[offsetInd]=yFine;
				offsetInd++;
			}
		}
		VLXindptrOut->data[offsetSep]=offsetInd;
		offsetSep++;
	}


	return 0;
}

int refine_signal(long int PartitionXAddr, TDoubleMatrix *signal, TDoubleMatrix *signalFine, int lTop) {
	THierarchicalPartition *partitionX=(THierarchicalPartition*) PartitionXAddr;
	partitionX->signal_refine_double(signal->data,signalFine->data,lTop);
	return 0;
}

int propagate_signalFunction(long int PartitionXAddr,
		TInteger64Matrix *signalPointer, int lTop, int lBottom, int mode) {

	THierarchicalPartition *partitionX=(THierarchicalPartition*) PartitionXAddr;
	double **signal=(double**) signalPointer->data;
	partitionX->signal_propagate_double(signal,lTop,lBottom,mode);

	return 0;
}

int propagate_costFunction(long int PartitionXAddr, long int PartitionYAddr,
		TInteger64Matrix *signalPointer, int lTop, int lBottom, int mode) {

	THierarchicalPartition *partitionX=(THierarchicalPartition*) PartitionXAddr;
	THierarchicalPartition *partitionY=(THierarchicalPartition*) PartitionYAddr;
	double **signal=(double**) signalPointer->data;

	THierarchicalProductSignal<double> *HierarchicalProduct=new THierarchicalProductSignal<double>(partitionX,partitionY);
	HierarchicalProduct->signal_propagate(signal,lTop,lBottom,mode);

	delete HierarchicalProduct;

	return 0;
}

int check_dualConstraints(long int PartitionXAddr, long int PartitionYAddr,
		TInteger64Matrix *cPointer, TInteger64Matrix *alphaPointer, TInteger64Matrix *betaPointer,
		int lTop, int lBottom, double slack,
		TInteger32Matrix *varListSpecs, TInteger64Matrix *varListPointer) {

	THierarchicalPartition *partitionX=(THierarchicalPartition*) PartitionXAddr;
	THierarchicalPartition *partitionY=(THierarchicalPartition*) PartitionYAddr;
	double **c=(double**) cPointer->data;
	double **alpha=(double**) alphaPointer->data;
	double **beta=(double**) betaPointer->data;

	THierarchicalProductSignal<double> *HierarchicalProduct=new THierarchicalProductSignal<double>
		(partitionX,partitionY);

	TVarListHandler *result;
	result=HierarchicalProduct->check_dualConstraints(c,alpha,beta,lTop,lBottom,slack);

	delete HierarchicalProduct;
	varListSpecs->data[0]=result->res;
	varListSpecs->data[1]=result->total;
	varListPointer->data[0]=(long int) result;
	return 0;
}

int check_dualConstraints_adaptive(long int PartitionXAddr, long int PartitionYAddr,
		TInteger64Matrix *cPointer, TInteger64Matrix *alphaPointer, TInteger64Matrix *betaPointer,
		int lTop, int lBottom, TInteger64Matrix *slackXPointer, TInteger64Matrix *slackYPointer,
		TInteger32Matrix *varListSpecs, TInteger64Matrix *varListPointer) {

	THierarchicalPartition *partitionX=(THierarchicalPartition*) PartitionXAddr;
	THierarchicalPartition *partitionY=(THierarchicalPartition*) PartitionYAddr;
	double **c=(double**) cPointer->data;
	double **alpha=(double**) alphaPointer->data;
	double **beta=(double**) betaPointer->data;
	double **slackX=(double**) slackXPointer->data;
	double **slackY=(double**) slackYPointer->data;

	THierarchicalProductSignal<double> *HierarchicalProduct=new THierarchicalProductSignal<double>
		(partitionX,partitionY);

	TVarListHandler *result;
	result=HierarchicalProduct->check_dualConstraints_adaptive(c,alpha,beta,lTop,lBottom,slackX,slackY);

	delete HierarchicalProduct;
	varListSpecs->data[0]=result->res;
	varListSpecs->data[1]=result->total;
	varListPointer->data[0]=(long int) result;
	return 0;
}


int Tools_Collect_VarList(TInteger32Matrix *indices, TInteger32Matrix *indptr,
		long int varListPointer, int doDelete) {
	TVarListHandler *varList;
	varList=(TVarListHandler*) varListPointer;
	varList->writeToCSRIndexList(indices->data,indptr->data);
	if (doDelete>0) {
		delete varList;
	}
	return 0;
}

int get_CEffectiveValues(TInteger32Matrix *indices, TInteger32Matrix *indptr,
		TDoubleMatrix *data,
		TDoubleMatrix *c, TDoubleMatrix *alpha, TDoubleMatrix *beta) {
	int x,y,xres,yres,yIndex;

	xres=c->dimensions[0];
	yres=c->dimensions[1];

	for(x=0;x<xres;x++) {
		for(yIndex=indptr->data[x];yIndex<(int) indptr->data[x+1];yIndex++) {
			y=indices->data[yIndex];
			data->data[yIndex]=c->data[x*yres+y]-alpha->data[x]-beta->data[y];
		}
	}
	return 0;
}

}
