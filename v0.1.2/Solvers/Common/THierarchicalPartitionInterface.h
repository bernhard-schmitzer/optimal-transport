#ifndef THierarchicalPartitionInterface_H_
#define THierarchicalPartitionInterface_H_

#include<stdlib.h>
#include<PythonTypes.h>
#include"THierarchicalPartition.h"

extern "C" {

int HierarchicalPartitionsCreate(int nLayers, int dim, TInteger64Matrix *PartitionPointer);
int HierarchicalPartitionsCreateLayer(long int PartitionAddr, int layerId,  int nCells, TInteger32Matrix *parent);
int HierarchicalPartitionsCreateLayerAddSublists(long int PartitionAddr, int layerId, int cellId,
	TInteger32Matrix *children, TInteger32Matrix *leaves);

int HierarchicalPartitionsGetSignalMass(long int PartitionAddr, TDoubleMatrix *mu, TInteger64Matrix *muLayers);

int HierarchicalPartitionsClose(long int PartitionAddr);

}

#endif /* THierarchicalPartitionInterface_H_ */
