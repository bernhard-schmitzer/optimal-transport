#ifndef TGeometry_H_
#define TGeometry_H_

#include<cstdio>
#include<Common/PythonTypes.h>
#include<Common/Tools.h>


class TGeometry_Euclidean {
	public:
	static double dist(const double * const posX, const double * const posY, const int dim) {
		double result=0;
		for(int i=0;i<dim;i++) {
			result+=std::pow(posX[i]-posY[i],2);
		}
		return std::sqrt(result);
	}
	static void geodesic(const double * const posX, const double * const posY,
			double * const posZ, const int dim, const double t) {

		for(int i=0;i<dim;i++) {
			posZ[i]=(1-t)*posX[i]+t*posY[i];
		}

	}
};

class TGeometry_Hyperbolic {
	public:
	const double scale;
	const double scaleSqr;
	
	const int bisections;
	
	TGeometry_Hyperbolic(const double _scale, const int _bisections) : scale(_scale), scaleSqr(_scale*_scale), bisections(_bisections) {};
	TGeometry_Hyperbolic(const double _scale) : TGeometry_Hyperbolic(_scale, 16) {};
	
	double dist(const double * const posX, const double * const posY, const int dim) const {
		double u2=EUCL_innerProduct(posX,posX,dim);
		double v2=EUCL_innerProduct(posY,posY,dim);
		double uv=EUCL_innerProduct(posX,posY,dim);
		double omega=(std::sqrt((scaleSqr+u2)*(scaleSqr+v2))-uv)/scaleSqr;
		if(omega<=1) {
			return 0.;
		}
		return scale*std::acosh(omega);
	}
	void geodesic(const double * const posX, const double * const posY,
			double * const posZ, const int dim, const double t) const {

		// for now only have parametrization of geodesic that is not arclength
		// function seems hard to inverse
		// for now sufficient: determine correct parameter via repeated bisection
		
		double u2=EUCL_innerProduct(posX,posX,dim);
		double v2=EUCL_innerProduct(posY,posY,dim);
		double uv=EUCL_innerProduct(posX,posY,dim);
		double omega=(std::sqrt((scaleSqr+u2)*(scaleSqr+v2))-uv)/scaleSqr;

		double sigma,tau;
		
		if(omega>1) {
			// determine appropriate parameter sigma via bisection
			
			double dist=std::acosh(omega); // full distance value, to compare against (ignore outer scaling here)
			
			sigma=0.5;
			for(int i=0;i<bisections;i++) {
				tau=std::sqrt(1+(std::pow(omega,2)-1)*std::pow(sigma,2))-omega*sigma;
				//printf("%d\t%e\t%e\n",i,sigma,tau);
				double omegaPath=(std::sqrt(
					(scaleSqr+u2)*
					(scaleSqr+(std::pow(sigma,2)*u2+std::pow(tau,2)*v2+2*tau*sigma*uv))
					)-(sigma*u2+tau*uv))/scaleSqr;
				double distPath=std::acosh(omegaPath);
				if(distPath>t*dist) {
					//printf("\t%f\t%f\t+\n",distPath/dist,t);
					sigma+=pow(0.5,i+2);
				} else{
					//printf("\t%f\t%f\t-\n",distPath/dist,t);
					sigma-=pow(0.5,i+2);
				}
			}
			tau=std::sqrt(1+(std::pow(omega,2)-1)*std::pow(sigma,2))-omega*sigma;
		} else {
			sigma=1.;
			tau=0.;
		}

		//printf("final result: %e\t%e\n",sigma,tau);
		for(int i=0;i<dim;i++) {
			posZ[i]=sigma*posX[i]+tau*posY[i];
		}

	}
};


#endif
