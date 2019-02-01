#ifndef Model_OT_H_
#define Model_OT_H_

#include<cstring>

#include<Common/PythonTypes.h>
#include<Common/Tools.h>

#include<Common/Models/TGeometry.h>

template<class TGeometry>
TParticleContainer ModelOT_Interpolate(const TSparsePosContainer& couplingData,
		const double * const posX, const double * const posY,
		const int dim, const double t, const TGeometry& geometry);



#endif
