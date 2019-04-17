#include"THierarchyBuilder.h"
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


THierarchyBuilder::THierarchyBuilder(double *_points, int _nPoints, int _dim,
		int _childMode, int partitionDepth) :
		points(_points), nPoints(_nPoints), dim(_dim), childMode(_childMode) {
	
	if (partitionDepth>0) {
		setBox();
		setRoot();
		// refine
		for(int i=0;i<partitionDepth-1;i++) {
			refine();
		}
	}
	addAtomicLayer();
	
	}

void THierarchyBuilder::setBox() {
	boxLo.resize(dim);
	boxHi.resize(dim);
	
	for(int i=0;i<dim;i++) {
		boxLo[i]=min(points,nPoints,dim,i)-boxTolerance;
		boxHi[i]=max(points,nPoints,dim,i)+boxTolerance;
	}
}
		

void THierarchyBuilder::setRoot() {
	// build top layer with only one node
	layers.resize(1);
	layers[0].nodes.resize(1);
	// aux variable for more convenient addressing of top node
	THierarchyBuilderNode *node;
	node=&(layers[0].nodes[0]);
	
	// set children to size zero
	node->children.resize(0);
	// set parent to 0
	node->parent=0;
	
	// set leaves to all points
	node->leaves.resize(nPoints);
	for(int i=0;i<nPoints;i++) {
		node->leaves[i]=i;
	}
	// set pos code to trivial pos code
	node->posCode.resize(dim);
	for(int i=0;i<dim;i++) {
		node->posCode[i]=0;
	}
	
}


void THierarchyBuilder::refine() {
	// create new finest layer
	layers.resize(layers.size()+1);
	int parentLayerId=layers.size()-2;
	
	if(childMode==CM_Grid) {
		// create empty layer of all grid located nodes
		int nChildNodes=pow(2,dim*(parentLayerId+1));
		layers[parentLayerId+1].nodes.resize(nChildNodes);
	}
	
	// aux vector to convert posCodes into global grid index
	std::vector<int> posVector(dim);
	for(int i=0;i<dim;i++) {
		posVector[i]=pow(2,(parentLayerId+1)*(dim-1-i));
	}
	
	for(int parentNodeId=0;parentNodeId<(int) layers[parentLayerId].nodes.size();parentNodeId++) {
		std::vector< std::vector<int> > childrenPosCodes(getChildrenPosCodes(parentLayerId,parentNodeId));
		std::vector< std::vector<int> > childrenLeaves(getChildrenLeaves(parentLayerId,parentNodeId,childrenPosCodes));
		// go over each children list
		for(int i=0;i<(int) childrenLeaves.size();i++) {
			int childNodeIndex;
			if(childMode==CM_Grid) {
				childNodeIndex=0;
				for(int k=0;k<dim;k++) {
					childNodeIndex+=childrenPosCodes[i][k]*posVector[k];
				}
			}
			if(childMode==CM_Tree) {
				// if leaves list is empty, skip this node. do not create a child
				if(childrenLeaves[i].size()==0) {
					continue;
				}
				// else: add empty node to layer
				layers[parentLayerId+1].nodes.resize(layers[parentLayerId+1].nodes.size()+1);
				childNodeIndex=layers[parentLayerId+1].nodes.size()-1;
			}
			// initialize corresponding child node
			THierarchyBuilderNode *node;
			node=&(layers[parentLayerId+1].nodes[childNodeIndex]);
			node->parent=parentNodeId;
			// pos code
			node->posCode=childrenPosCodes[i];
			// leaves list
			node->leaves=childrenLeaves[i];
			// add new child node to child list of parent
			layers[parentLayerId].nodes[parentNodeId].children.push_back(childNodeIndex);
		}
	}
}

std::vector<std::vector<int> > THierarchyBuilder::getChildrenPosCodes(int layerId, int nodeId) {
	// for all children of a node, construct the corresponding refined posCodes.

	// pointer to parent node for convenience
	THierarchyBuilderNode *node=&(layers[layerId].nodes[nodeId]);
	// total number of potential children
	int nSublists=pow(2,dim);
	
	// create empty children posCode lists
	std::vector<std::vector<int> > result;
	result.resize(nSublists);

	// iterate over all possible child posCodes
	
	for(int i=0;i<nSublists;i++) {
		// allocate childs posCode list
		result[i].resize(dim);
		
		// set pos code
		// get relativ pos code for child
		getRelPosCodeFromIndex(i,dim,&result[i][0]);
		
		// get absolute pos code
		getOffsetPosCode(&result[i][0],&node->posCode[0],dim);						
	}
	
	return result;
}

std::vector<std::vector<int> > THierarchyBuilder::getChildrenLeaves(int layerId, int nodeId,
		const std::vector<std::vector<int> >& childrenPosCodes) {
	// for all children of a given node, construct the corresponding leaves lists

	// pointer to parent node for convenience
	THierarchyBuilderNode *node=&(layers[layerId].nodes[nodeId]);
	// number of leaves of parent
	int nLeaves=node->leaves.size();
	// total number of potential children
	int nSublists=pow(2,dim);

	// create available list for all leaves. mark all leaves available at beginning.
	std::vector<bool> available(nLeaves);
	for(int i=0;i<nLeaves;i++) {
		available[i]=true;
	}
	
	// create empty children leave lists
	std::vector<std::vector<int> > result;
	result.resize(nSublists);

	// iterate over all possible child posCodes
	for(int i=0;i<nSublists;i++) {
		// go over all leaves and check if they are in sub-box
		for(int j=0;j<nLeaves;j++) {
			if(available[j]) {
				int leaveId=node->leaves[j];
				if(isInBox(points+leaveId*dim,&childrenPosCodes[i][0],dim,layerId+1)) {
					available[j]=false;
					result[i].push_back(leaveId);
				}
			}
		}
	}
	
	return result;
}

void THierarchyBuilder::getRelPosCodeFromIndex(int index, int dim, int *posCode) {
	for(int j=0;j<dim;j++) {
		posCode[j]=(index % ( (int) pow(2,dim-j) ) ) / pow(2,dim-j-1);
	}
}

void THierarchyBuilder::getOffsetPosCode(int *relPosCode, int *parentPosCode, int dim) {
	for(int i=0;i<dim;i++) {
		relPosCode[i]+=parentPosCode[i]*2;
	}
}

bool THierarchyBuilder::isInBox(double *coord, const int * const posCode, int dim, int layerId) {
	// go through each dim and check upper and lower bound of box
	for(int i=0;i<dim;i++) {
		double width=(boxHi[i]-boxLo[i])*pow(2,-layerId); // width of box in that dimension
		if((coord[i]<boxLo[i]+width*posCode[i]-boxTolerance) || (coord[i]>boxLo[i]+width*(posCode[i]+1)+boxTolerance)) {
			return false;
		}
	}
	return true;
}

void THierarchyBuilder::addAtomicLayer() {
	// create new finest layer
	layers.resize(layers.size()+1);
	int layerId=layers.size()-1;
	
	layers[layerId].nodes.resize(nPoints);
	
	if(layerId>0) {
		// write parent values
		// iterate over coarse nodes
		for(int i=0;i<(int) layers[layerId-1].nodes.size();i++) {
			// iterate over all leaves, set children to leaves and set parent value at all children
			layers[layerId-1].nodes[i].children=layers[layerId-1].nodes[i].leaves;
			for(int j=0;j<(int) layers[layerId-1].nodes[i].leaves.size();j++) {
				layers[layerId].nodes[
					layers[layerId-1].nodes[i].leaves[j]
					].parent=i;
			}
		}
	}
}


double THierarchyBuilder::max(double *x, int n, int step, int offset) {
	double result=x[offset];
	for(int i=1;i<n;i++) {
		if(result<x[i*step+offset]) {
			result=x[i*step+offset];
		}
	}
	return result;
}

double THierarchyBuilder::min(double *x, int n, int step, int offset) {
	double result=x[offset];
	for(int i=1;i<n;i++) {
		if(result>x[i*step+offset]) {
			result=x[i*step+offset];
		}
	}
	return result;
}



THierarchicalPartition* THierarchyBuilder::convert() {
	THierarchicalPartition *result;
	result=new THierarchicalPartition(layers.size(),dim);
	
	// construct layers
	for(int i=0;i<(int) layers.size();i++) {
		// construct layer with empty nodes
		result->layers[i]=new TPartitionLayer();
		result->layers[i]->initializeEmpty(layers[i].nodes.size());
		
		// export parent data
		result->layers[i]->parent=(int*) malloc(sizeof(int)*layers[i].nodes.size());
		for(int j=0;j<(int) layers[i].nodes.size();j++) {
			result->layers[i]->parent[j]=layers[i].nodes[j].parent;
		}
		
		// if not on finest level: export children and leaves data
		if(i<(int) layers.size()-1) {
			for(int j=0;j<(int) layers[i].nodes.size();j++) {
				// leaves
				result->layers[i]->leaves[j]=(int*) malloc(sizeof(int)*layers[i].nodes[j].leaves.size());
				result->layers[i]->nLeaves[j]=layers[i].nodes[j].leaves.size();
				for(int k=0;k<(int) layers[i].nodes[j].leaves.size();k++) {
					result->layers[i]->leaves[j][k]=layers[i].nodes[j].leaves[k];
				}
				// children
				result->layers[i]->children[j]=(int*) malloc(sizeof(int)*layers[i].nodes[j].children.size());
				result->layers[i]->nChildren[j]=layers[i].nodes[j].children.size();
				for(int k=0;k<(int) layers[i].nodes[j].children.size();k++) {
					result->layers[i]->children[j][k]=layers[i].nodes[j].children[k];
				}
			}				
		}
	}
	
	return result;
}

double** THierarchyBuilder::allocateDoubleSignal(int sigdim, int nLayers) {
	// allocates a double signal with sigdim entries per point, for nLayers layers (starting at 0)
	double **result;
	int _nLayers;
	if(nLayers==0) {
		_nLayers=layers.size();
	} else {
		_nLayers=nLayers;
	}
	result=(double**) malloc(sizeof(double*)*layers.size());
	for(int i=0;i< _nLayers;i++) {
		result[i]=(double*) malloc(sizeof(double)*layers[i].nodes.size()*sigdim);
	}
	return result;
}

void THierarchyBuilder::freeSignal(double **signal, int nLayers) {
	for(int i=0;i<nLayers;i++) {
		free(signal[i]);
	}
	free(signal);
}

void THierarchyBuilder::getSignalPos(double **signal) {
	int i;
	// non-atomic layers
	for(i=0;i<(int) layers.size()-1;i++) {
		for(int k=0;k<dim;k++) {
			double width=(boxHi[k]-boxLo[k])*pow(2,-i); // width of box in that dimension
			for(int j=0;j<(int) layers[i].nodes.size();j++) {
				signal[i][j*dim+k]=boxLo[k]+width*(layers[i].nodes[j].posCode[k]+0.5);
			}
		}
	}
	// atomic layer
	i=layers.size()-1;
	for(int j=0;j<nPoints*dim;j++) {
		signal[i][j]=points[j];
	}
}

int* THierarchyBuilder::getDimH(int *finestDims) {
	int *result=(int*) malloc(sizeof(int)*layers.size()*dim);
	
	int i;
	// non-atomic layers
	for(i=0;i< (int) layers.size()-1;i++) {
		for(int j=0;j<dim;j++) {
			result[i*dim+j]=pow(2,i);
		}
	}
	
	// atomic layer
	i=layers.size()-1;
	for(int j=0;j<dim;j++) {
		result[i*dim+j]=finestDims[j];
	}
	return result;
	
}

int* THierarchyBuilder::getResH() {
	int *result=(int*) malloc(sizeof(int)*layers.size());
	for(int i=0;i<(int) layers.size();i++) {
		result[i]=(int) layers[i].nodes.size();
	}
	return result;

}


double** THierarchyBuilder::getSignalRadii() {
	return getSignalRadiiAxis(NULL,0);
}

double** THierarchyBuilder::getSignalRadiiAxis(int *axis, int naxis) {
	// only compute radii for certain subset of axis.
	// axis is int array that indicates which axis are meant
	// naxis gives length
	
	// for convenience: when calling with (*axis=NULL,naxis=0) all axis are used.

	// allocate memory for radii
	double **result=allocateDoubleSignal(1, layers.size()-1);
	
	// compute diagonal of full bounding box
	double diagonal=0;
	
	if(naxis==0) {
		// sum squared lenghts along all axis
		for(int i=0;i<dim;i++) {
			diagonal+=pow(boxHi[i]-boxLo[i],2);
		}
	} else{
		// sum squared lenghts along selected axis
		for(int i=0;i<naxis;i++) {
			diagonal+=pow(boxHi[axis[i]]-boxLo[axis[i]],2);
		}
	}
	
	// take sqrt, factor 0.5: turn diameter into radius
	diagonal=0.5*sqrt(diagonal);
	
	for(int iLayer=0; iLayer< (int) layers.size()-1;iLayer++) {
		// scale down according to layer number: 0.5 per layer
		double subdiagonal=diagonal*pow(0.5,iLayer);
		// write value to all signal entries at given layer
		for(int iNode=0; iNode< (int) layers[iLayer].nodes.size();iNode++) {
			result[iLayer][iNode]=subdiagonal;
		}
	}
	
	return result;
}

