#include"tools.h"

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

