#include"Tools.h"

void doubleArrayCopy(double *a, double *b, int n) {
	int i;
	for(i=0;i<n;i++) {
		b[i]=a[i];
	}
}

void doubleArrayScale(double *a, double b, int n) {
	int i;
	for(i=0;i<n;i++) {
		a[i]=a[i]*b;
	 }
}


double doubleArrayMin(double *data, int res) {
	// compute minimum of array
	double result=data[0];
	for(int i=1;i<res;i++) {
		if(result>data[i]) {
			result=data[i];
		}
	}
	return result;
}

double doubleArrayMax(double *data, int res) {
	// compute maximum of array
	double result=data[0];
	for(int i=1;i<res;i++) {
		if(result<data[i]) {
			result=data[i];
		}
	}
	return result;
}


int intArrayMin(int *data, int res) {
	// compute minimum of array
	int result=data[0];
	for(int i=1;i<res;i++) {
		if(result>data[i]) {
			result=data[i];
		}
	}
	return result;
}

int intArrayMax(int *data, int res) {
	// compute maximum of array
	int result=data[0];
	for(int i=1;i<res;i++) {
		if(result<data[i]) {
			result=data[i];
		}
	}
	return result;
}


void freeTDoubleMatrix(TDoubleMatrix *a) {
	free(a->dimensions);
	free(a->data);
	free(a);
}


////////////////////////////////

// a few geometric aux functions that are often used
double EUCL_innerProduct(const double * const a, const double * const b, const int n) {
	// computes inner product a.b, dimension given by n
	double result=0;
	int i;
	for(i=0;i<n;i++) {
		result+=a[i]*b[i];
	}
	return result;
}

void EUCL_lincomb(const double * const a, const double * const b, double * const c, const double sa, const double sb, const int n) {
	// stores sa*a+sb*b to c, dimension given by n
	int i;
	for(i=0;i<n;i++) {
		c[i]=sa*a[i]+sb*b[i];
	}
}

double EUCL_lincombSqr(const double * const a, const double * const b, const double sa, const double sb, const int n) {
	// returns |sa*a+sb*b|^2, dimension given by n
	int i;
	double result,buffer;
	result=0;
	for(i=0;i<n;i++) {
		buffer=sa*a[i]+sb*b[i];
		result+=buffer*buffer;
	}
	return result;
}

double EUCL_len(const double * const a, const int n) {
	int i;
	double result=0;
	for(i=0;i<n;i++) {
		result+=a[i]*a[i];
	}
	return pow(result,0.5);
}


///////////////////////////////////////////////////////
// read raw data from file

template <typename V>
std::vector<V> readFile(const char* filename) {
	std::ifstream input( filename, std::ios::binary );
	// copies all data into buffer
	std::vector<char> buffer((
			std::istreambuf_iterator<char>(input)), 
			(std::istreambuf_iterator<char>()));

	int len=(int) buffer.size()/sizeof(V);
	std::vector<V> result(len);
	std::memcpy( result.data(), buffer.data(), buffer.size() );
	return result;
}


template <typename V>
int writeFile(const char* filename, const V * const data, const int len) {
	std::ofstream output( filename, std::ios::binary );
	char* dataByte=(char*) data;
	int lenByte=len*sizeof(V);
	output.write(dataByte,lenByte);
	
	return 0;
}


template std::vector<double> readFile<double>(const char* filename);
template int writeFile<double>(const char* filename, const double* const data, const int len);

