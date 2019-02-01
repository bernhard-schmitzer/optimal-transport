#ifndef Model_WFR_H_
#define Model_WFR_H_

#include<algorithm>

#include<Common/PythonTypes.h>
#include<Common/Tools.h>

#include<Common/Models/OT.h>
#include<Common/Models/TGeometry.h>

template<class TGeometry>
TParticleContainer ModelWFR_Interpolate(const TSparsePosContainer& couplingData,
		const double * const muXEff, const double * const muYEff,
		const double * const muX, const double * const muY,
		const double * const posX, const double * const posY,
		const int dim, const double t, const double kappa,
		const TGeometry& geometry);


#endif
