#include"GridTools.h"


void GridToolsGetStrides(int dim, int *dims, int *strides) {
	int d;
	strides[dim-1]=1;
	for(d=dim-2;d>=0;d--) {
		strides[d]=strides[d+1]*dims[d+1];
	}
}

int GridToolsGetIdFromPos(int dim, int *pos, int *strides) {
	int d, result;
	result=0;
	for(d=0;d<dim;d++) {
		result+=pos[d]*strides[d];
	}
	return result;
}

void GridToolsGetPosFromId(int dim, int id, int *pos, int *strides) {
	int d;
	pos[0]=id/strides[0];
	for(d=1;d<dim;d++) {
		pos[d]=(id%strides[d-1])/strides[d];
	}
}

//////////////////////////////

void GridToolsGetNeighbours(int dim, int *dims, TVarListHandler *neighbours) {
		GridToolsGetNeighbours_Torus(dim,dims,0,neighbours);
}

void GridToolsGetNeighbours_Torus(int dim, int *dims, int torusDim, TVarListHandler *neighbours) {
	// dim: dimension of grid
	// dims: points per axis
	// torusDim: dimensions 0 to torusDim-1 are considered cyclical.
		// dimensions torusDim to dim-1 are considered standard Euclidean
	
	int *strides;
	strides=(int*) malloc(sizeof(int)*dim);
	int *pos;
	pos=(int*) malloc(sizeof(int)*dim);

	GridToolsGetStrides(dim,dims,strides);

	GridToolsGetNeighbours_Torus_iterateXVariables(neighbours, pos, dims, strides, dim, torusDim, 0);

	free(strides);
	free(pos);
}

void GridToolsGetNeighbours_Torus_iterateXVariables(TVarListHandler *neighbours, int *pos, int *dims, int *strides,
		int dim, int torusDim, int d) {
	int x;
	if(d<dim) {
		// as long as we have not reached the last axis, launch iteration along next axis
		for(x=0;x<dims[d];x++) {
			pos[d]=x;
			GridToolsGetNeighbours_Torus_iterateXVariables(neighbours,pos,dims,strides,dim,torusDim,d+1);
		}
	} else {
		// once a specific grid point is reached, add points
		//addVariables_Shields(xVars,xMap,xPos);
		int xId,d;
		xId=GridToolsGetIdFromPos(dim,pos,strides);
		for(d=0;d<dim;d++) {
			if(pos[d]>0) {
				neighbours->addToLine(xId,xId-strides[d]);
			} else {
				if(d<torusDim) {
					neighbours->addToLine(xId,xId+strides[d]*(dims[d]-1));
				}
			}
			if(pos[d]<dims[d]-1) {
				neighbours->addToLine(xId,xId+strides[d]);
			} else {
				if(d<torusDim) {
					neighbours->addToLine(xId,xId-strides[d]*(dims[d]-1));
				}
			}
		}

	}
}




int GridToolsGetTotalPoints(int depth, int *dimensions) {
	int totalSize=1;
	for(int d=0;d<depth;d++) {
		totalSize*=dimensions[d];
	}
	return totalSize;
}

double* GridToolsGetGrid(int depth, int *dimensions) {
	double *result;
	// total number of grid points
	int totalSize=GridToolsGetTotalPoints(depth,dimensions);
	
	// allocate grid points
	result=(double*) malloc(sizeof(double)*totalSize*depth);
	
	// set grid points
	// iterate along each dimension. set coordinate of that dimension for all matching points
	for(int d=0;d<depth;d++) {
		// set combined dimension cardinality variables
		// for all coarser dimensions
		int cardCoarse=GridToolsGetTotalPoints(d,dimensions);

		// for finer dimensions
		int cardFine=GridToolsGetTotalPoints(depth-d-1,dimensions+d+1);
		
		for(int i=0;i<cardCoarse;i++) {
			for(int j=0;j<dimensions[d];j++) {
				for(int k=0;k<cardFine;k++) {
					result[i*dimensions[d]*cardFine*depth+j*cardFine*depth+k*depth+d]=j;
				}
			}
		}
	}
	return result;
}

TDoubleMatrix* GridToolsGetGridMatrix(int depth, int *dimensions) {
	TDoubleMatrix *result=(TDoubleMatrix*) malloc(sizeof(TDoubleMatrix));
	result->data=GridToolsGetGrid(depth,dimensions);
	result->depth=2;
	result->dimensions=(int*) malloc(sizeof(int)*2);
	result->dimensions[0]=GridToolsGetTotalPoints(depth,dimensions);
	result->dimensions[1]=depth;
	return result;
	
}


void GridToolsShift(double *pos, double *shift, double scale, int nPoints, int dim) {
	for(int i=0;i<nPoints;i++) {
		for(int j=0;j<dim;j++) {
			pos[i*dim+j]+=scale*shift[j];
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////



double MeasureToolsTruncateMeasure(double *muX, int xres, double measureScale) {
	// truncate a single measure to integer values, , by dividing through measureScale and rounding.
	// return total truncated mass
	double result=0;
	for(int i=0;i<xres;i++) {
		muX[i]=round(muX[i]/measureScale);
		result+=muX[i];
	}
	return result;
}

int MeasureToolsTruncateMeasures(double *muX, double *muY, int xres, int yres, double measureScale) {
	// truncate two measures to integer values, by dividing through measureScale and rounding.
	// small differences in total mass are then compensated by adjusting measure with lower mass
	double muXsum,muYsum;
	muXsum=MeasureToolsTruncateMeasure(muX,xres,measureScale);
	muYsum=MeasureToolsTruncateMeasure(muY,yres,measureScale);
	
	// compensate mass difference by iterating over measure with lower mass and increasing entries until balance is restored
	double *muLess;
	int resLess;
	int increments;
	if(muXsum<muYsum) {
		muLess=muX;
		resLess=xres;
		increments=round(muYsum-muXsum);
	} else {
		muLess=muY;
		resLess=yres;
		increments=round(muXsum-muYsum);
	}
	int offset=0;
	while(increments>0) {
		muLess[offset]+=1.;
		increments--;
		offset++;
		if(offset>=resLess) {
			offset=0;
		}
	}

	// if somee entries have been rounded to zero, throw error
	if(doubleArrayMin(muX,xres)<=0.) {
		return ERR_PREP_TRUNC_MUXNEG;
	}
	if(doubleArrayMin(muY,yres)<=0.) {
		return ERR_PREP_TRUNC_MUYNEG;
	}
	
	// now multiply both measures once more by measureScale
	for(int x=0;x<xres;x++) {
		muX[x]*=measureScale;
	}
	for(int y=0;y<yres;y++) {
		muY[y]*=measureScale;
	}

	// return 0 to indicate successful truncation
	return 0;	
}


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


int MeasureToolsTruncateMeasureToInt(double *muXIn, int *muXOut, int xres, double measureScale) {
	// truncate a single measure to integer values, by dividing through measureScale and rounding.
	// return total truncated mass
	int result=0;
	for(int i=0;i<xres;i++) {
		muXOut[i]=lround(muXIn[i]/measureScale);
		result+=muXOut[i];
	}
	return result;
}

int MeasureToolsTruncateMeasuresToInt(
		double *muXIn, double *muYIn, int *muXOut, int *muYOut,
		int xres, int yres, double measureScale) {
	
	// truncate two measures to integer values, by dividing through measureScale and rounding.
	// small differences in total mass are then compensated by adjusting measure with lower mass
	int muXsum,muYsum;
	muXsum=MeasureToolsTruncateMeasureToInt(muXIn,muXOut,xres,measureScale);
	muYsum=MeasureToolsTruncateMeasureToInt(muYIn,muYOut,yres,measureScale);
	
	// compensate mass difference by iterating over measure with lower mass and increasing entries until balance is restored
	int *muLess;
	int resLess;
	int increments;
	if(muXsum<muYsum) {
		muLess=muXOut;
		resLess=xres;
		increments=muYsum-muXsum;
	} else {
		muLess=muYOut;
		resLess=yres;
		increments=muXsum-muYsum;
	}
	int offset=0;
	while(increments>0) {
		muLess[offset]+=1;
		increments--;
		offset++;
		if(offset>=resLess) {
			offset=0;
		}
	}

	// if somee entries have been rounded to zero, throw error
	if(intArrayMin(muXOut,xres)<=0) {
		return ERR_PREP_TRUNC_MUXNEG;
	}
	if(intArrayMin(muYOut,yres)<=0) {
		return ERR_PREP_TRUNC_MUYNEG;
	}

	// return 0 to indicate successful truncation
	return 0;	
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


template<int dim>
int ProjectInterpolation(const double * const pos, const double * const mass,
		double * const data,
		const int nParticles, const int* const res) {
	// rasterizes list of particles with mass to underlying grid
	// use simple bi-linear interpolation between four nearest grid points
	// pos: array of positions of particles (assumed to be array of 2d coordinates)
	// mass: array of particle masses
	// data: array of output image, is assumed to be initialized. new data will simply be added
	// nParticles: number of particles
	// res: array of grid dimensions
	
	double weight[dim];
	int gridPos[dim];
	
	
	for(int i=0;i<nParticles;i++) {
		// compute grid positions and mass distribution ratios for rasterization
		for(int d=0;d<dim;d++) {
			gridPos[d]=std::max(0,std::min(res[d]-1,(int) pos[dim*i+d]));
			if((pos[dim*i+d]>=0) && (pos[dim*i+d]<=res[d]-1)) {
				weight[d]=1-std::abs(((double) gridPos[d])-pos[dim*i+d]);
			} else {
				weight[d]=1;
			}
		}

		// add mass to appropriate grid positions
		switch(dim) {
			case 2:
				ProjectInterpolation_SingleParticle2D(mass[i], data, res, gridPos, weight);
				break;

			default:
				return ERR_BASE_NOTIMPLEMENTED;			
		}
	}
	
	return 0;
	
}


template<int dim>
int ProjectInterpolationVector(const double * const pos, const double * const mass,
		double * const data, double * const dataMass,
		const int nParticles, const int* const res) {
	// rasterizes list of particles with mass to underlying grid
	// vector version: only first two coordinates of pos are spatial, rest is part of signal that is to be added linearly, weighed by mass
	
	// use simple bi-linear interpolation between four nearest grid points
	// pos: array of positions of particles (assumed to be array of 2d coordinates)
	// mass: array of particle masses
	// data: array of output image, IS ASSUMED TO BE ZERO INITIALLY (!!!!).
	// dataMass: mass of output image, IS ASSUMED TO BE ZERO INITIALLY (!!!!)
	// nParticles: number of particles
	// res: array of grid dimensions, last dimension is sigDim.
	
	// fullDim: full dimension of pos coordinates (i.e. 2 spatial dimensions + signal dimension)
	
	
	double weight[dim];
	int gridPos[dim];
	const int sigDim=res[dim]; // last axis of res stores sigDim
	const int fullDim=dim+sigDim;
	
	
	for(int i=0;i<nParticles;i++) {
		// compute grid positions and mass distribution ratios for rasterization
		for(int d=0;d<dim;d++) {
			gridPos[d]=std::max(0,std::min(res[d]-1,(int) pos[fullDim*i+d]));
			if((pos[fullDim*i+d]>=0) && (pos[fullDim*i+d]<=res[d]-1)) {
				weight[d]=1-std::abs(((double) gridPos[d])-pos[fullDim*i+d]);
			} else {
				weight[d]=1;
			}
		}

		// add mass and signal to appropriate grid positions
		switch(dim) {
			case 2:
				// rasterize mass
				ProjectInterpolation_SingleParticle2D(mass[i], dataMass, res, gridPos, weight);
				// rasterize signal
				ProjectInterpolation_SingleParticle2D_Vector(
						mass[i], pos+(fullDim*i+dim) /* pointer to signal of that particle */,
						data, res,
						gridPos, weight);
				
				break;

			default:
				return ERR_BASE_NOTIMPLEMENTED;			
		}
	}
	
	// final normalization of signal: divide signal at pixel by total mass at pixel
	
	// compute nr of total pixels, by multiplying spatial components of res
	int spatRes=1;
	for(int d=0;d<dim;d++) {
		spatRes*=res[d];
	}
	
	// now loop over pixels
	for(int i=0;i<spatRes;i++) {
		if(dataMass[i]>0) {
			// loop over signal dimensions
			for(int s=0;s<sigDim;s++) {
				data[i*sigDim+s]=data[i*sigDim+s]/dataMass[i];
			}
		}
	}
	
	return 0;
	
}


void ProjectInterpolation_SingleParticle2D(
		const double mass,
		double * const data, const int * const res,
		const int * const gridPos, const double * const weight) {

	data[gridPos[0]*res[1]+gridPos[1]]+=weight[0]*weight[1]*mass;
	if(gridPos[0]<res[0]-1) {
		data[(gridPos[0]+1)*res[1]+gridPos[1]]+=(1-weight[0])*weight[1]*mass;
	}
	if(gridPos[1]<res[1]-1) {
		data[gridPos[0]*res[1]+(gridPos[1]+1)]+=weight[0]*(1-weight[1])*mass;
	}
	if((gridPos[0]<res[0]-1) && (gridPos[1]<res[1]-1)) {
		data[(gridPos[0]+1)*res[1]+(gridPos[1]+1)]+=(1-weight[0])*(1-weight[1])*mass;
	}

}

void ProjectInterpolation_SingleParticle2D_Vector(
		const double mass, const double * const sig,
		double * const data, const int * const res,
		const int * const gridPos, const double * const weight) {

	const int dim=2; // number of spatial dimensions
	const int sigDim=res[dim]; // last axis of res stores sigDim

	// iterate over signal dimensions
	for(int s=0;s<sigDim;s++) {
		data[gridPos[0]*res[1]*sigDim+gridPos[1]*sigDim+s]+=weight[0]*weight[1]*mass*sig[s];
		if(gridPos[0]<res[0]-1) {
			data[(gridPos[0]+1)*res[1]*sigDim+gridPos[1]*sigDim+s]+=
				(1-weight[0])*weight[1]*mass*sig[s];
		}
		if(gridPos[1]<res[1]-1) {
			data[gridPos[0]*res[1]*sigDim+(gridPos[1]+1)*sigDim+s]+=
				weight[0]*(1-weight[1])*mass*sig[s];
		}
		if((gridPos[0]<res[0]-1) && (gridPos[1]<res[1]-1)) {
			data[(gridPos[0]+1)*res[1]*sigDim+(gridPos[1]+1)*sigDim+s]+=
				(1-weight[0])*(1-weight[1])*mass*sig[s];
		}
	}

	
}


template int ProjectInterpolation<2>(const double * const pos, const double * const mass,
		double * const data,
		const int nParticles, const int* const res);

template int ProjectInterpolationVector<2>(const double * const pos, const double * const mass,
		double * const data, double * const dataMass,
		const int nParticles, const int* const res);



