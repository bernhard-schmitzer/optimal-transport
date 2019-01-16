#include"THierarchicalPartitionInterface.h"

int HierarchicalPartitionsCreate(int nLayers, int dim, TInteger64Matrix *PartitionPointer) {
	THierarchicalPartition *partition;

	partition=new THierarchicalPartition(nLayers,dim);

	PartitionPointer->data[0]=(long int) partition;

	return 0;
}

int HierarchicalPartitionsCreateLayer(long int PartitionAddr, int layerId,  int nCells,
	TInteger32Matrix *parent) {

	THierarchicalPartition *partition=(THierarchicalPartition*) PartitionAddr;

	partition->layers[layerId]=new TPartitionLayer();
	partition->layers[layerId]->initializeEmpty(nCells);
	partition->layers[layerId]->parent=parent->data;
	
	return 0;
}

int HierarchicalPartitionsCreateLayerAddSublists(long int PartitionAddr, int layerId, int cellId,
	TInteger32Matrix *children, TInteger32Matrix *leaves) {

	THierarchicalPartition *partition=(THierarchicalPartition*) PartitionAddr;

	partition->layers[layerId]->children[cellId]=children->data;
	partition->layers[layerId]->leaves[cellId]=leaves->data;
	partition->layers[layerId]->nChildren[cellId]=children->dimensions[0];
	partition->layers[layerId]->nLeaves[cellId]=leaves->dimensions[0];

	return 0;	
}

int HierarchicalPartitionsClose(long int PartitionAddr) {
	THierarchicalPartition *partition=(THierarchicalPartition*) PartitionAddr;
	delete partition;
	return 0;
}


int HierarchicalPartitionsGetSignalMass(long int PartitionAddr, TDoubleMatrix *mu, TInteger64Matrix *muLayers) {
	THierarchicalPartition *partition=(THierarchicalPartition*) PartitionAddr;
	partition->computeHierarchicalMasses(mu->data,(double**) muLayers->data);
	return 0;
}

