#include"Geometry_Sphere.h"

	
int TMultiScaleSetupSingleSphere::SetupProjectPoints() {
	for(int layer=0;layer<nLayers;layer++) {
		SetupProjectPoints_Array(posH[layer], HP->layers[layer]->nCells, dim);
	}
	return 0;
}

int TMultiScaleSetupSingleSphere::SetupProjectPoints_Array(double *pos, int n, int dim) {
	double len;
	for(int i=0;i<n;i++) {
		len=EUCL_len(pos+dim*i,dim);
		if(len<sphereCenterTolerance) {
			// set to north pole
			pos[dim*i]=1.;
			for(int j=1;j<dim;j++) {
				pos[dim*i+j]=0.;
			}
		} else {
			// rescale
			for(int j=0;j<dim;j++) {
				pos[dim*i+j]=pos[dim*i+j]/len;
			}
		}
	}
	return 0;
}

double TMultiScaleSetupSingleSphere::SphereDistance(double *posx, double *posy) {
	double result;
	result=EUCL_innerProduct(posx, posy, dim);
	if(result>=1.) { return 0.; }
	if(result<=-1.) { return M_PI; }
	result=std::acos(result);
	return result;
}

int TMultiScaleSetupSingleSphere::SetupRadii() {
	radii=HP->signal_allocate_double(0,HP->nLayers-2);
		
	// iterate from fine to coarse
	// on coarsest layer: radii are zero, no data allocated for this
	// on each other layer: compute maximum over (distance+child radius) over all children
	for(int layer=HP->nLayers-2;layer>=0;layer--) {
		
		int nCells=HP->layers[layer]->nCells;
		for(int x=0;x<nCells;x++) {
			int nChildren=HP->layers[layer]->nChildren[x];
			double currentMax=0, currentRadius;
			for(int yId=0;yId<nChildren;yId++) {
				int y=HP->layers[layer]->children[x][yId];
				currentRadius=SphereDistance(posH[layer]+x*dim,posH[layer+1]+y*dim);
				if(layer<HP->nLayers-2) {
					currentRadius+=radii[layer+1][y];	
				}
				currentMax=std::max(currentMax,currentRadius);
			}
			radii[layer][x]=currentMax;
		}
	}
	
	return 0;
}




///////////////////////////////////////////////////////////////////////////////////////////////////	
// HierarchicalCostFunctionProvider
///////////////////////////////////////////////////////////////////////////////////////////////////	


double THierarchicalCostFunctionProvider_Sphere::getCostAsym(int layerX, int x, int layerY, int y) {
	double result;
	
	// compute sphere distance from inner product
	result=EUCL_innerProduct(xPos[layerX]+(x*posDim), yPos[layerY]+(y*posDim), posDim);
	if(result>=1.) {
		result=0.;
	} else {
		if(result<=-1.) {
			result=M_PI;
		} else {
			result=std::acos(result);
		}
	}	

	
	// if not at finest layer, need to compute lower bound
	if(layerX<layerBottom) {
		result-=xRadii[layerX][x];
	}
	if(layerY<layerBottom) {
		result-=yRadii[layerY][y];
	}
	if(result<0) {
		result=0;
	}

	return std::pow(result,p);
}


