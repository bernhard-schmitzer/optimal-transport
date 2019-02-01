#ifndef PythonTypes_H_
#define PythonTypes_H_


struct TDoubleMatrix {
	double *data;
	int depth;
	int *dimensions;
	};

struct TInteger32Matrix {
	int *data;
	int depth;
	int *dimensions;
	};

struct TInteger64Matrix {
	long int *data;
	int depth;
	int *dimensions;
	};

struct TCSRDoubleMatrix {
	double *data;
	int *indices;
	int *indptr;
	int xres,yres;
	int total;
	};


struct TCSRInteger32Matrix {
	int *data;
	int *indices;
	int *indptr;
	int xres,yres;
	int total;
	};


#endif
