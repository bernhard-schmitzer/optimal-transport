#include"Geometry_Sphere.h"

	
int TMultiScaleSetupSphere::SetupProjectPoints() {
	for(int layer=0;layer<nLayers;layer++) {
		SetupProjectPoints_Array(posXH[layer], HPX->layers[layer]->nCells, dim);
		SetupProjectPoints_Array(posYH[layer], HPY->layers[layer]->nCells, dim);
	}
	return 0;
}

int TMultiScaleSetupSphere::SetupProjectPoints_Array(double *pos, int n, int dim) {
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

double TMultiScaleSetupSphere::SphereDistance(double *posx, double *posy) {
	double result;
	result=EUCL_innerProduct(posx, posy, dim);
	if(result>=1.) { return 0.; }
	if(result<=-1.) { return M_PI; }
	result=std::acos(result);
	return result;
}

int TMultiScaleSetupSphere::SetupRadii_Single(THierarchicalPartition *HPA, double **aPosH, double ***aRadii) {
	*aRadii=HPA->signal_allocate_double(0,HPA->nLayers-2);
		
	// iterate from fine to coarse
	// on coarsest layer: radii are zero, no data allocated for this
	// on each other layer: compute maximum over (distance+child radius) over all children
	for(int layer=HPA->nLayers-2;layer>=0;layer--) {
		
		int nCells=HPA->layers[layer]->nCells;
		for(int x=0;x<nCells;x++) {
			int nChildren=HPA->layers[layer]->nChildren[x];
			double currentMax=0, currentRadius;
			for(int yId=0;yId<nChildren;yId++) {
				int y=HPA->layers[layer]->children[x][yId];
				currentRadius=SphereDistance(aPosH[layer]+x*dim,aPosH[layer+1]+y*dim);
				if(layer<HPA->nLayers-2) {
					currentRadius+=(*aRadii)[layer+1][y];	
				}
				currentMax=std::max(currentMax,currentRadius);
			}
			(*aRadii)[layer][x]=currentMax;
		}
	}
	
	return 0;
}

int TMultiScaleSetupSphere::SetupRadii() {
	SetupRadii_Single(HPX,posXH,&xRadii);
	SetupRadii_Single(HPY,posYH,&yRadii);
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


