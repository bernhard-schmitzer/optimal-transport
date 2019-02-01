#include"OT.h"


template<class TGeometry>
TParticleContainer ModelOT_Interpolate(const TSparsePosContainer& couplingData,
		const double * const posX, const double * const posY,
		const int dim, const double t, const TGeometry& geometry) {
	// compute displacement interpolation in geometry specified by TGeometry
	
	// couplingData contains coupling data as TSparsePosContainer.
	//  couplingData.mass: particle masses,
	//  couplingData.posStart and posEnd: indices of initial and final position
	// posX, posY: spatial position of marginal particles
	// dim: spatial dimension
	// t: time in [0,1] for which to compute interpolation snapshot
	
	// returns TParticleContainer
	//  mass: masses of particles
	//  pos: spatial positions of particles
	
	TParticleContainer result;
	// same number of particles
	result.nParticles=couplingData.nParticles;
	// same masses, simply copy
	result.mass=couplingData.mass;

	// allocate spatial position array
	result.pos.resize(result.nParticles*dim);
	// iterate over particles and dimensions and set result to linear interpolation
	for(int i=0;i<result.nParticles;i++) {
		geometry.geodesic(posX+(couplingData.posStart[i]*dim), posY+(couplingData.posEnd[i]*dim),result.pos.data()+(i*dim),dim,t);
	}
	
	return result;
}


template TParticleContainer ModelOT_Interpolate<TGeometry_Euclidean>(const TSparsePosContainer& couplingData,
		const double * const posX, const double * const posY,
		const int dim, const double t, const TGeometry_Euclidean& geometry);

template TParticleContainer ModelOT_Interpolate<TGeometry_Hyperbolic>(const TSparsePosContainer& couplingData,
		const double * const posX, const double * const posY,
		const int dim, const double t, const TGeometry_Hyperbolic& geometry);


