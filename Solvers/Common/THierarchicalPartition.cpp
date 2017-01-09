#include"THierarchicalPartition.h"

TPartitionLayer::TPartitionLayer() {
	nCells=0;
	parent=NULL;
	children=NULL;
	leaves=NULL;
	nChildren=NULL;
	nLeaves=NULL;
}

TPartitionLayer::~TPartitionLayer() {
	if (children!=NULL) {
		free(children);
	}
	if (leaves!=NULL) {
		free(leaves);
	}
	if (nChildren!=NULL) {
		free(nChildren);
	}
	if (nLeaves!=NULL) {
		free(nLeaves);
	}
}

void TPartitionLayer::initializeEmpty(int _nCells) {
	nCells=_nCells;
	// allocate space for children and leave lists
	
	children=(int**) malloc(sizeof(int*)*nCells);
	leaves=(int**) malloc(sizeof(int*)*nCells);

	nChildren=(int*) malloc(sizeof(int)*nCells);
	nLeaves=(int*) malloc(sizeof(int)*nCells);

	for(int i=0;i<nCells;i++) {
		children[i]=NULL;
		leaves[i]=NULL;
		nChildren[i]=0;
		nLeaves[i]=0;
	}
}


THierarchicalPartition::THierarchicalPartition(int _nLayers, int _dim) {
	nLayers=_nLayers;
	dim=_dim;
	
	layers = (TPartitionLayer**) malloc(sizeof(TPartitionLayer*)*nLayers);
	for (int i=0;i<nLayers;i++) {
		layers[i]=NULL;
	}
}

THierarchicalPartition::~THierarchicalPartition() {
	for(int i=0;i<nLayers;i++) {
		if(layers[i]!=NULL) {
			delete layers[i];
		}
	}
	free(layers);
}


void THierarchicalPartition::computeHierarchicalMasses(double *mu, double **muLayers) {
	int l,x,y;
	// copy lowest layer
	for(x=0;x<layers[nLayers-1]->nCells;x++) {
		muLayers[nLayers-1][x]=mu[x];
	}
	// now go up through higher layers, each time summing mass of children
	l=nLayers-2;
	while(l>=0) {
		// go through all cells of current layer
		for(x=0;x<layers[l]->nCells;x++) {
			// reset mass of current cell
			muLayers[l][x]=0;
			// now go through children of this cell
			for(y=0;y<layers[l]->nChildren[x];y++) {
				muLayers[l][x]+=muLayers[l+1][
					layers[l]->children[x][y]
					];
			}
		}
		l--;
	}
}

double** THierarchicalPartition::signal_allocate_double(int lTop, int lBottom) {
	double **result;
	result=(double**) malloc(sizeof(double*)*(lBottom-lTop+1));
	for(int i=0;i<lBottom-lTop+1;i++) {
		result[i]=(double*) malloc(sizeof(double)* layers[lTop+i]->nCells);
	}
	
	return result;
}

void THierarchicalPartition::signal_free_double(double **signal, int lTop, int lBottom) {
	for(int i=0;i<=lBottom-lTop;i++) {
		free(signal[i]);
	}
	free(signal);

}

void THierarchicalPartition::signal_propagate_double(double **signal, int lTop, int lBottom, int mode) {
	double newValue,value;
	// iterate through all involved layers (except the finest one), from bottom up
	for(int i=lBottom-1;i>=lTop;i--) {
		// go through all cells on current layer
		for(int j=0;j<layers[i]->nCells;j++) {
			// initialize signal with value of first child
			value=signal[i-lTop+1][layers[i]->children[j][0] ];
			// iterate over rest of children
			for(int k=1;k<layers[i]->nChildren[j];k++) {
				newValue=signal[i-lTop+1][layers[i]->children[j][k] ];
				if( ((mode==MODE_MAX) && (newValue>value)) || ((mode==MODE_MIN) && (newValue<value))) {
					value=newValue;
				}
			}
			signal[i-lTop][j]=value;
		}
	}
}

void THierarchicalPartition::signal_refine_double(double *signal, double *signalFine, int lTop) {
	int x,xFineIndex,xFine;
	// iterate over top layer
	for(x=0;x<layers[lTop]->nCells;x++) {
		// iterate over children
		for(xFineIndex=0;xFineIndex<layers[lTop]->nChildren[x];xFineIndex++) {
			xFine=layers[lTop]->children[x][xFineIndex];
			// set signal to parent value
			signalFine[xFine]=signal[x];
		}
	}
}

template<typename T>
THierarchicalProductSignal<T>::THierarchicalProductSignal(THierarchicalPartition *_partitionX, THierarchicalPartition *_partitionY) {
	partitionX=_partitionX;
	partitionY=_partitionY;
}


template<typename T>
void THierarchicalProductSignal<T>::signal_propagate(T **signal, int lTop, int lBottom, int mode) {
	T newValue,value;
	TPartitionLayer *layerX, *layerY;
	int yres;
	int x,y,xFineIndex,yFineIndex,xFine,yFine;
	// iterate through all involved layers (except the finest one), from bottom up
	for(int i=lBottom-1;i>=lTop;i--) {
		layerX=partitionX->layers[i];
		layerY=partitionY->layers[i];
		yres=partitionY->layers[i+1]->nCells;
		// go through all cells on current layer in X
		for(x=0;x<(int) layerX->nCells;x++) {
			// go through all cells on current layer in Y
			for(y=0;y<(int) layerY->nCells;y++) {
				value=0;
				// iterate over children
				for(xFineIndex=0;xFineIndex<(int) layerX->nChildren[x];xFineIndex++) {
					xFine=layerX->children[x][xFineIndex];
					for(yFineIndex=0;yFineIndex<(int) layerY->nChildren[y];yFineIndex++) {
						yFine=layerY->children[y][yFineIndex];
						newValue=signal[i-lTop+1][xFine*yres + yFine ];
						if(((xFineIndex==0) && (yFineIndex==0))
								|| ((mode==MODE_MAX) && (newValue>value))
								|| ((mode==MODE_MIN) && (newValue<value))
								) {
							value=newValue;
						}
					}
				}
				signal[i-lTop][x*layerY->nCells+y]=value;
			}
		}
	}
}

template<typename T>
TVarListHandler* THierarchicalProductSignal<T>::check_dualConstraints(T **_c, T **_alpha, T **_beta, int lTop, int lBottom, T _slack) {
	TVarListHandler* result;
	int xres,x,yres,y;
	
	// store params to temporariy global variables
	c=_c;
	alpha=_alpha;
	beta=_beta;
	slack=_slack;
	
	// create varList
	xres=partitionX->layers[lBottom]->nCells;
	result=new TVarListHandler();
	result->setupEmpty(xres);
	
	// set global varList
	varList=result;
	
	// initialize iterations
	xres=partitionX->layers[lTop]->nCells;
	yres=partitionY->layers[lTop]->nCells;
	for(x=0;x<xres;x++) {
		for(y=0;y<yres;y++) {
			check_dualConstraints_iterateTile(lTop, x, y, lBottom);
		}
	}
	
	// reset global variables
	c=NULL;
	alpha=NULL;
	beta=NULL;
	varList=NULL;
	
	// return
	return result;
}

template<typename T>
void THierarchicalProductSignal<T>::check_dualConstraints_iterateTile(int l, int x, int y, int lBottom) {
	// analyze cell (x,y) on layer l. That is go through all its children on layer l+1, check constraints there.
	// if violated either refine (if l+1<lBottom) or add variable
	int yres,xFine,yFine,xFineIndex,yFineIndex;
	
	//xres=partitionX->layers[l+1]->nCells;
	yres=partitionY->layers[l+1]->nCells;
	
	for(xFineIndex=0;xFineIndex<(int) partitionX->layers[l]->nChildren[x];xFineIndex++) {
		xFine=partitionX->layers[l]->children[x][xFineIndex];
		
		for(yFineIndex=0;yFineIndex<(int) partitionY->layers[l]->nChildren[y];yFineIndex++) {
			yFine=partitionY->layers[l]->children[y][yFineIndex];
			if(c[l+1][xFine*yres+yFine]-alpha[l+1][xFine]-beta[l+1][yFine]<=slack) {
				if(l+1==lBottom) {
					varList->addToLine(xFine,yFine,false);
				} else{
					check_dualConstraints_iterateTile(l+1, xFine, yFine, lBottom);
				}
			}
		}
	}
	
}



template<typename T>
TVarListHandler* THierarchicalProductSignal<T>::check_dualConstraints_adaptive(T **_c, T **_alpha, T **_beta, int lTop, int lBottom,
		T **_slackOffsetX, T **_slackOffsetY) {
	TVarListHandler* result;
	int xres,x,yres,y;
	
	// store params to temporariy global variables
	c=_c;
	alpha=_alpha;
	beta=_beta;
	slackOffsetX=_slackOffsetX;
	slackOffsetY=_slackOffsetY;
	
	// create varList
	xres=partitionX->layers[lBottom]->nCells;
	result=new TVarListHandler();
	result->setupEmpty(xres);
	
	// set global varList
	varList=result;
	
	// initialize iterations
	xres=partitionX->layers[lTop]->nCells;
	yres=partitionY->layers[lTop]->nCells;
	for(x=0;x<xres;x++) {
		for(y=0;y<yres;y++) {
			check_dualConstraints_adaptive_iterateTile(lTop, x, y, lBottom);
		}
	}
	
	// reset global variables
	c=NULL;
	alpha=NULL;
	beta=NULL;
	varList=NULL;
	slackOffsetX=NULL;
	slackOffsetY=NULL;
	
	// return
	return result;
}

template<typename T>
void THierarchicalProductSignal<T>::check_dualConstraints_adaptive_iterateTile(int l, int x, int y, int lBottom) {
	// analyze cell (x,y) on layer l. That is go through all its children on layer l+1, check constraints there.
	// if violated either refine (if l+1<lBottom) or add variable
	int yres,xFine,yFine,xFineIndex,yFineIndex;
	double slackValue;
	
	//xres=partitionX->layers[l+1]->nCells;
	yres=partitionY->layers[l+1]->nCells;
	
	for(xFineIndex=0;xFineIndex<(int) partitionX->layers[l]->nChildren[x];xFineIndex++) {
		xFine=partitionX->layers[l]->children[x][xFineIndex];
		
		for(yFineIndex=0;yFineIndex<(int) partitionY->layers[l]->nChildren[y];yFineIndex++) {
			yFine=partitionY->layers[l]->children[y][yFineIndex];
			slackValue=c[l+1][xFine*yres+yFine]-alpha[l+1][xFine]-beta[l+1][yFine];
			if((slackValue<=slackOffsetX[l+1][xFine]) || (slackValue<=slackOffsetY[l+1][yFine])) {
				//cout << l+1 << "\t" << xFine << "\t" << yFine << "\t" << slackValue << "\t"
				//		<< slackOffsetX[l+1][xFine] << "\t" << slackOffsetY[l+1][yFine] << "\n" << endl;
				if(l+1==lBottom) {
					varList->addToLine(xFine,yFine,false);
				} else{
					check_dualConstraints_adaptive_iterateTile(l+1, xFine, yFine, lBottom);
				}
			}
		}
	}
	
}


template class THierarchicalProductSignal<double>;

