#include"TShieldGenerator_Models.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TShieldGeneratorTree_Sphere
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<class TReferenceTree>
TShieldGeneratorTree_SpherePrototype<TReferenceTree>::TShieldGeneratorTree_SpherePrototype(int _dim,
		THierarchicalPartition* _yPartition, double** _yPos, double** _yRadii,
		int _lBottom, int _lTop,
		int *_xresH, double **_xposH, TVarListHandler **_xNeighboursH,
		double _p) :
				TReferenceTree(_dim, _yPartition, _yPos, _yRadii,
						_lBottom, _lTop,
						_xresH, _xposH, _xNeighboursH) {

	p=_p;
}



template<class TReferenceTree>
TShieldGeneratorTree_SpherePrototype<TReferenceTree>::~TShieldGeneratorTree_SpherePrototype() {
}


template<class TReferenceTree>
bool TShieldGeneratorTree_SpherePrototype<TReferenceTree>::checkConditionPlane(int xA, int x, int l,
		int yB, int y) {
	// xA: index of initial x
	// yB: index of final y
	// x,y: indices of shielding pair
	// l: hierarchy level. the higher, the finer the resolution.

	// inner products
	double IPxAy, IPxy, IPxAx, IPxAyB, IPxyB;
	// final shielding bound estimate
	double dist;

	IPxAy=EUCL_innerProduct(xpos+(xA*dim),yPos[lBottom]+(y*dim),dim);
	IPxy=EUCL_innerProduct(xpos+(x*dim),yPos[lBottom]+(y*dim),dim);


	IPxAx=EUCL_innerProduct(xpos+(xA*dim),xpos+(x*dim),dim);
	IPxAyB=EUCL_innerProduct(xpos+(xA*dim),yPos[l]+(yB*dim),dim);
	IPxyB=EUCL_innerProduct(xpos+(x*dim),yPos[l]+(yB*dim),dim);


	// compute r.h.s. of shielding condition
	dist=-(getPhi(acos(IPxAy))-getPhi(acos(IPxy)));

	if(l>=lBottom) {
		// if at finest level,
		dist+=getPhi(acos(IPxAyB))-getPhi(acos(IPxyB));
		return (dist>=shieldingTolerance);
	} else {

		double theta, phi;
		//double thetaS;
		double deltaTheta, deltaPhi;
		double subgrad;
		double thetaMin, phiMax;
		double deltaD;

		double cosDeltaTheta,sinThetaS,sinTheta;

		theta=acos(IPxAyB);
		//thetaS=acos(IPxAx);
		sinTheta=std::sqrt(1-IPxAyB*IPxAyB);
		sinThetaS=std::sqrt(1-IPxAx*IPxAx);

		deltaTheta=yRadii[l][yB];
		cosDeltaTheta=cos(deltaTheta);

		phi=acos(std::min(1.,std::max(-1.,(IPxyB-IPxAx*IPxAyB)/(sinTheta*sinThetaS ) ) ) );

		// compute delta phi


		double cosDeltaThetaSqr,cosThetaSqr;
		// if point lies close to one of the poles
		cosDeltaThetaSqr=cosDeltaTheta*cosDeltaTheta;
		cosThetaSqr=IPxAyB*IPxAyB;
		if(cosThetaSqr>=cosDeltaThetaSqr) {
			deltaPhi=M_PI;
		} else {
			deltaPhi=acos(sqrt((cosDeltaThetaSqr-cosThetaSqr)/(1-cosThetaSqr)));

		}


		thetaMin=std::max(0.,theta-deltaTheta);
		phiMax=std::min(M_PI,phi+deltaPhi);


		double d1,d2;
		d1=thetaMin;
		d2=acos(sinThetaS*sin(thetaMin)*cos(phiMax)+IPxAx*cos(thetaMin));

		/*
		// fixed square variant
		deltaD=d1*d1-d2*d2;
		dist+=deltaD;
		return (dist>shieldingTolerance);
		*/


		// general subgradient variant

		deltaD=(d1-d2);

		// get subgradient bound
		if(deltaD>0) {
			subgrad=getSubgrad(std::max(0.,acos(IPxyB)-deltaTheta));
		} else {
			subgrad=getSubgrad(std::min(M_PI,acos(IPxyB)+deltaTheta));
		}

		dist+=subgrad*deltaD;

		return (dist>shieldingTolerance);

	}
}


template class TShieldGeneratorTree_SpherePrototype<TShieldGeneratorTreeBase>;
//template class TShieldGeneratorTree_SpherePrototype<TShieldGeneratorTreeBase_Benchmark>;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TShieldGeneratorTree_Reflector
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<class TReferenceTree>
TShieldGeneratorTree_ReflectorPrototype<TReferenceTree>::TShieldGeneratorTree_ReflectorPrototype(int _dim,
		THierarchicalPartition* _yPartition, double** _yPos, double** _yRadii,
		int _lBottom, int _lTop,
		int *_xresH, double **_xposH, TVarListHandler **_xNeighboursH
		) :
				TReferenceTree(_dim, _yPartition, _yPos, _yRadii,
						_lBottom, _lTop,
						_xresH, _xposH, _xNeighboursH) {

}


template<class TReferenceTree>
TShieldGeneratorTree_ReflectorPrototype<TReferenceTree>::~TShieldGeneratorTree_ReflectorPrototype() {
}


template<class TReferenceTree>
bool TShieldGeneratorTree_ReflectorPrototype<TReferenceTree>::checkConditionPlane(int xA, int x, int l,
		int yB, int y) {
	// xA: index of initial x
	// yB: index of final y
	// x,y: indices of shielding pair
	// l: hierarchy level. the higher, the finer the resolution.

	// inner products
	double IPxAy, IPxy, IPxAx, IPxAyB, IPxyB;
	// final shielding bound estimate
	double dist;

	IPxAy=EUCL_innerProduct(xpos+(xA*dim),yPos[lBottom]+(y*dim),dim);
	IPxy=EUCL_innerProduct(xpos+(x*dim),yPos[lBottom]+(y*dim),dim);


	IPxAx=EUCL_innerProduct(xpos+(xA*dim),xpos+(x*dim),dim);
	IPxAyB=EUCL_innerProduct(xpos+(xA*dim),yPos[l]+(yB*dim),dim);
	IPxyB=EUCL_innerProduct(xpos+(x*dim),yPos[l]+(yB*dim),dim);


	// compute r.h.s. of shielding condition
	dist=-(getPhi(acos(IPxAy))-getPhi(acos(IPxy)));

	if(l>=lBottom) {
		// if at finest level,
		dist+=getPhi(acos(IPxAyB))-getPhi(acos(IPxyB));
		return (dist>=shieldingTolerance);
	} else {

		double theta, phi;
		//double thetaS;
		double deltaTheta, deltaPhi;
		double subgrad;
		double thetaMax, phiMin;
		double deltaD;

		double cosDeltaTheta,sinThetaS,sinTheta;

		theta=acos(IPxAyB);
		//thetaS=acos(IPxAx);
		sinTheta=std::sqrt(1-IPxAyB*IPxAyB);
		sinThetaS=std::sqrt(1-IPxAx*IPxAx);


		deltaTheta=yRadii[l][yB];
		cosDeltaTheta=cos(deltaTheta);

		phi=acos(std::min(1.,std::max(-1.,(IPxyB-IPxAx*IPxAyB)/(sinTheta*sinThetaS ) ) ) );

		// compute delte phi

		double cosDeltaThetaSqr,cosThetaSqr;
		// if point lies close to one of the poles
		cosDeltaThetaSqr=cosDeltaTheta*cosDeltaTheta;
		cosThetaSqr=IPxAyB*IPxAyB;
		if(cosThetaSqr>=cosDeltaThetaSqr) {
			deltaPhi=M_PI;
		} else {
			deltaPhi=acos(sqrt((cosDeltaThetaSqr-cosThetaSqr)/(1-cosThetaSqr)));

		}


		thetaMax=std::min(M_PI,theta+deltaTheta);
		phiMin=std::max(0.,phi-deltaPhi);


		double d1,d2;
		d1=thetaMax;
		d2=acos(sinThetaS*sin(thetaMax)*cos(phiMin)+IPxAx*cos(thetaMax));

		// general subgradient variant

		deltaD=(d1-d2);

		// get subgradient bound
		if(deltaD>0) {
			subgrad=getSubgrad(std::max(0.,acos(IPxyB)-deltaTheta));
		} else {
			subgrad=getSubgrad(std::min(M_PI,acos(IPxyB)+deltaTheta));
		}

		dist+=subgrad*deltaD;

		return (dist>shieldingTolerance);

	}
}

template class TShieldGeneratorTree_ReflectorPrototype<TShieldGeneratorTreeBase>;

