#ifndef GridTools_H_
#define GridTools_H_

#include<cmath>

#include<Common/TVarListHandler.h>
#include<Common/ErrorCodes.h>
#include<Common/PythonTypes.h>
#include<Common/Tools.h>

void GridToolsGetStrides(int dim, int *dims, int *strides);
int GridToolsGetIdFromPos(int dim, int *pos, int *strides);
void GridToolsGetPosFromId(int dim, int id, int *pos, int *strides);

void GridToolsGetNeighbours(int dim, int *dims, TVarListHandler *neighbours);
void GridToolsGetNeighbours_Torus(int dim, int *dims, int torusDim, TVarListHandler *neighbours);
void GridToolsGetNeighbours_Torus_iterateXVariables(TVarListHandler *neighbours, int *pos, int *dims, int *strides,
		int dim, int torusDim, int d);


int GridToolsGetTotalPoints(int depth, int *dimensions);
double* GridToolsGetGrid(int depth, int *dimensions);
TDoubleMatrix* GridToolsGetGridMatrix(int depth, int *dimensions);

void GridToolsShift(double *pos, double *shift, double scale, int nPoints, int dim);

double MeasureToolsTruncateMeasure(double *muX, int xres, double measureScale);
int MeasureToolsTruncateMeasures(double *muX, double *muY, int xres, int yres, double measureScale);
int MeasureToolsTruncateMeasureToInt(double *muXIn, int *muXOut, int xres, double measureScale);
int MeasureToolsTruncateMeasuresToInt(
		double *muXIn, double *muYIn, int *muXOut, int *muYOut,
		int xres, int yres, double measureScale);



template<int dim>
int ProjectInterpolation(const double * const pos, const double * const mass,
		double * const data,
		const int nParticles, const int* const res);

template<int dim>
int ProjectInterpolationVector(const double * const pos, const double * const mass,
		double * const data, double * const dataMass,
		const int nParticles, const int* const res);

void ProjectInterpolation_SingleParticle2D(
		const double mass,
		double * const data, const int * const res,
		const int * const gridPos, const double * const weight);

void ProjectInterpolation_SingleParticle2D_Vector(
		const double mass, const double * const sig,
		double * const data, const int * const res,
		const int * const gridPos, const double * const weight);

#endif

