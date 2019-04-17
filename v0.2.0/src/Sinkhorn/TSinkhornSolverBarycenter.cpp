#include"TSinkhornSolverBarycenter.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// standard Wasserstein barycenter
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TSinkhornSolverBarycenter::TSinkhornSolverBarycenter(
		int _nLayers,
		int *_nEpsList,
		double **_epsLists,
		int _layerCoarsest, int _layerFinest,
		TSinkhornSolverParameters _cfg,
		int _nMarginals,
		double *_weights,
		THierarchicalPartition **_HP, THierarchicalPartition *_HPZ,
		double ***_muH, double **_muZH,
		double ***_alphaH, double ***_betaH,		
		THierarchicalCostFunctionProvider **_costProvider
		) : TSinkhornSolverBase(_nLayers, _nEpsList, _epsLists, _layerCoarsest, _layerFinest, _cfg) {
	
	nMarginals=_nMarginals;
	weights=_weights;
	HP=_HP;
	HPZ=_HPZ;
	muH=_muH;
	muZH=_muZH;
	alphaH=_alphaH;
	betaH=_betaH;
	costProvider=_costProvider;
	
	kernelGenerator=NULL;

}

TSinkhornSolverBarycenter::~TSinkhornSolverBarycenter() {
	if(kernelGenerator!=NULL) {
		for(int i=0;i<nMarginals;i++) {
			delete kernelGenerator[i];
		}
		free(kernelGenerator);
	}
}

int TSinkhornSolverBarycenter::initialize() {
	// call basic init method
	TSinkhornSolverBase::initialize();
	
	kernel.resize(nMarginals);
	kernelT.resize(nMarginals);
	
	res.resize(nMarginals,0);
	mu.resize(nMarginals);
	alpha.resize(nMarginals);
	beta.resize(nMarginals);
	u.resize(nMarginals);
	v.resize(nMarginals);

	
	// initialize kernel generators
	kernelGenerator=(TSinkhornKernelGenerator**) malloc(sizeof(TSinkhornKernelGenerator*)*nMarginals);
	for(int i=0;i<nMarginals;i++) {
		kernelGenerator[i]=new TSinkhornKernelGenerator(
			HP[i], HPZ,
			muH[i][layer], muZH[layer],
			costProvider[i],
			layer
			);
	}

	return 0;
}


int TSinkhornSolverBarycenter::refineDuals(const int newLayer) {

	for(int i=0;i<nMarginals;i++) {		
		// refine dual variables
		if(newLayer>0) {
			HP[i]->signal_refine_double(alphaH[i][newLayer-1],alphaH[i][newLayer],newLayer-1);
			HPZ->signal_refine_double(betaH[i][newLayer-1],betaH[i][newLayer],newLayer-1);
		}
	}
	kernelValid=false;
	return 0;


}

int TSinkhornSolverBarycenter::changeLayer(const int newLayer) {
	TSinkhornSolverBase::changeLayer(newLayer);

	// things that only need to be done once (not once per marginal)	
	zres=HPZ->layers[layer]->nCells;
	muZ=muZH[layer];
		
	
	for(int i=0;i<nMarginals;i++) {		
		// adjust cost function provider
		costProvider[i]->setLayerBottom(layer);
	
		// adjust kernel generator
		kernelGenerator[i]->layerBottom=layer;
		kernelGenerator[i]->muX=muH[i][layer];
		kernelGenerator[i]->muY=muZH[layer];
		
		// set short cut variables for current layer
		res[i]=HP[i]->layers[layer]->nCells;
		mu[i]=muH[i][layer];
		alpha[i]=alphaH[i][layer];
		beta[i]=betaH[i][layer];

		// setup marginal variables
		u[i]=TMarginalVector::Constant(res[i],1.);
		v[i]=TMarginalVector::Constant(zres,1.);
	
	}
	
	return 0;

}


int TSinkhornSolverBarycenter::checkAbsorb(const double maxValue) {
	for(int i=0;i<nMarginals;i++) {
		if((u[i].maxCoeff()>maxValue) || (v[i].maxCoeff()>maxValue)) {
			return MSG_ABSORB_REITERATE;
		}
	}
	return 0;
}

int TSinkhornSolverBarycenter::absorb() {
	for(int i=0;i<nMarginals;i++) {
		SinkhornAbsorbScaling(HP[i], alphaH[i], u[i], layer, eps);
		SinkhornAbsorbScaling(HPZ, betaH[i], v[i], layer, eps);
	}
	return 0;	
}

int TSinkhornSolverBarycenter::generateKernel() {
	for(int i=0;i<nMarginals;i++) {
		// set slack for effective cost based on thresh for kernel entry
		kernel[i]=kernelGenerator[i]->generate(eps,-eps*log(cfg.truncation_thresh),false,false);

		kernelT[i]=kernel[i].transpose();
		eprintf("\t\tkernel entries: %ld\n",kernel[i].nonZeros());
	}
	return 0;
}


int TSinkhornSolverBarycenter::refineKernel() {
	for(int i=0;i<nMarginals;i++) {
		kernel[i]=kernelGenerator[i]->refine(kernel[i],eps);
		kernelT[i]=kernel[i].transpose();
		eprintf("\t\tkernel entries: (refinement) %ld\n",kernel[i].nonZeros());	
	}
	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Standard optimal transport model

int TSinkhornSolverBarycenter::iterate(const int n) {
	// standard sinkhorn iterations
	
	// aux variables
	std::vector<TMarginalVector> logConvList(nMarginals);
	
	for(int i=0;i<n;i++) { // iterate over number of iterations
	
		// u-update
		for(int j=0;j<nMarginals;j++) {
			u[j]=Eigen::Map<TMarginalVector>(mu[j],res[j]).cwiseQuotient(kernel[j]*v[j]);
			if(!u[j].allFinite()) {
				eprintf("\tNAN in u[%d]\n",j);
				return MSG_NANSCALING;
			}
		}



		// v-update	
		
		// first compute averages (that handle the "communication" between all marginals)
		TMarginalVector logConvAvg=TMarginalVector::Constant(zres,0.);
		TMarginalVector absorbedAvg=TMarginalVector::Constant(zres,0.);
		TMarginalVector convolution;

		// indicator of where one of the kernel columns is empty and thus the iteration is not defined
		//TMarginalVectorInt zeroSet=TMarginalVectorInt::Constant(zres,0);
		TMarginalArrayInt zeroSet=TMarginalArrayInt::Constant(zres,0);
		
		for(int j=0;j<nMarginals;j++) { // iterate over number of marginals
			convolution=kernelT[j]*u[j]; // convolution with (transposed) kernel
			
			// add entries where convolution is zero (i.e. kernel col is empty) to zero set
			zeroSet+=(convolution.array()<DBL_ZEROCOLTOLERANCE).cast<int>();
			
			logConvList[j]=convolution.array().log().matrix();
			logConvAvg+=weights[j]*logConvList[j];
			absorbedAvg+=weights[j]*Eigen::Map<TMarginalVector>(beta[j],zres);


		}		


		// apply changes to scaling factors
		for(int j=0;j<nMarginals;j++) { // iterate over number of marginals
		
			// compute naive new values of v[j]
			TMarginalVector vNew=(logConvAvg-absorbedAvg/eps-logConvList[j]).array().exp().matrix();
			// only apply new value outside of zero set. on zero set keep old values
			v[j]=(zeroSet>0).select(v[j],vNew);
			if(!v[j].allFinite()) {
				eprintf("\tNAN in v[%d]\n",j);
				return MSG_NANSCALING;
			}
		}
		
	}	
		
		
	return 0;
}


int TSinkhornSolverBarycenter::getError(double * const result) {
	// return sum of L1 differences of barycenter marginals from their mean
	std::vector<TMarginalVector> innerMarginals(nMarginals);
	TMarginalVector meanMarginal=TMarginalVector::Constant(zres,0.);
	for(int i=0;i<nMarginals;i++) {
		// compute 2nd marginal of i-th transport plan (for i-th reference measure)
		innerMarginals[i]=v[i].cwiseProduct(kernelT[i]*u[i]);
		// add to meanMarginal
		meanMarginal+=innerMarginals[i];
	}
	// normalize mean
	meanMarginal=meanMarginal/nMarginals;
	(*result)=0;
	for(int i=0;i<nMarginals;i++) {
		(*result)+=(meanMarginal-innerMarginals[i]).lpNorm<1>();
	}

	return 0;

}

//double TSinkhornSolverBarycenter::scoreTransportCost() {
//	// evaluates linear transport cost term
//	double result=0;
//	
//	// iterate over elements of kernel
//	for (int x=0; x<kernel.outerSize(); x++) {
//		for (TKernelMatrix::InnerIterator it(kernel,x); it; ++it) {
//			int y=it.col();
//			// cost contribution: cost * pi(x,y), where pi(x,y) = u(x)*v(y)*kernel(x,y) (the latter is stored in iterator value)
//			result+=costProvider->getCost(layer, x, y)*u[x]*v[y]*it.value();
//		}
//	}	
//	return result;
//}

//double TSinkhornSolverBarycenter::scorePrimalUnreg() {
//	// evalutes unregularized primal functional.
//	// ignore marginal constraints, simply return transport cost term
//	return scoreTransportCost();
//}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// KL marginal version

TSinkhornSolverBarycenterKLMarginals::TSinkhornSolverBarycenterKLMarginals(
		int _nLayers,
		int *_nEpsList,
		double **_epsLists,
		int _layerCoarsest, int _layerFinest,
		TSinkhornSolverParameters _cfg,
		int _nMarginals,
		double *_weights,
		THierarchicalPartition **_HP, THierarchicalPartition *_HPZ,
		double ***_muH, double **_muZH,
		double ***_alphaH, double ***_betaH,		
		THierarchicalCostFunctionProvider **_costProvider,
		double _kappa
		) : TSinkhornSolverBarycenter(
			_nLayers, _nEpsList, _epsLists,
			_layerCoarsest, _layerFinest,
			_cfg,
			_nMarginals,
			_weights,
			_HP, _HPZ, _muH, _muZH, _alphaH, _betaH,		
			_costProvider)
		{	
	
	kappa=_kappa;	
}	


int TSinkhornSolverBarycenterKLMarginals::iterate(const int n) {
	// aux variable that stores convolutions with kernel
	// code is interspersed with commented code from the python implementation
	
	std::vector<TMarginalArray> convolution(nMarginals);
	for(int i=0;i<n;i++) { // iterate over number of iterations

		// u-update
		for(int j=0;j<nMarginals;j++) {
			convolution[j]=(kernel[j]*v[j]).array();
			TMarginalArrayInt zeroSet=(convolution[j]<DBL_ZEROCOLTOLERANCE).cast<int>();
			TMarginalVector uNew=(
					(Eigen::Map<TMarginalArray>(mu[j],res[j])/convolution[j]).pow(kappa/(kappa+eps))
					*(Eigen::Map<TMarginalArray>(alpha[j],res[j])/(-kappa-eps)).exp()
					).matrix();
			u[j]=(zeroSet>0).select(u[j],uNew);
			if(!u[j].allFinite()) {
				eprintf("\tNAN in u[%d]\n",j);
				return MSG_NANSCALING;
			}
		}

	
		// v-update	
		
		// phi is aux variable that contains some form of average over all convolutions
		TMarginalArray phi=TMarginalArray::Constant(zres,0.);

		
		for(int j=0;j<nMarginals;j++) { // iterate over number of marginals
			convolution[j]=(kernelT[j]*u[j]).array(); // convolution with (transposed) kernel
			//cout << "conv" << j << endl << convolution[j] << endl << endl;
			
			phi+=weights[j]*convolution[j].pow(eps/(eps+kappa))
					*(Eigen::Map<TMarginalArray>(beta[j],zres)/(-kappa-eps)).exp();
			//phi+=weights[i]*((a[i])**(eps/(prefac+eps)))*np.exp(-alphaList[n+i]/(prefac+eps))

		}
		//phi=phi**(prefac/eps)
		phi=phi.pow(kappa/eps);
		//cout << "phi" << endl << phi << endl << endl;
			


		// apply changes to scaling factors
		for(int j=0;j<nMarginals;j++) { // iterate over number of marginals
		
			// compute naive new values of v[j]
			TMarginalVector vNew=(
					phi/(convolution[j].pow(kappa/(kappa+eps)))
					*(Eigen::Map<TMarginalArray>(beta[j],zres)/(-kappa-eps)).exp()
					).matrix();
			//scalingList[n+i][:]=phi/(a[i]**(prefac/(prefac+eps)))*np.exp(-alphaList[n+i]/(prefac+eps))
					
			// only apply new value outside of zero set. on zero set keep old values
			TMarginalArrayInt zeroSet=(convolution[j]<DBL_ZEROCOLTOLERANCE).cast<int>();

			v[j]=(zeroSet>0).select(v[j],vNew);
			if(!v[j].allFinite()) {
				eprintf("\tNAN in v[%d]\n",j);				
				return MSG_NANSCALING;
			}
			v[j]=v[j].cwiseMax(DBL_MINSCALING);
		}
		
	}	
		
		
	return 0;
}

int TSinkhornSolverBarycenterKLMarginals::getError(double * const result) {
	double scorePrimal,scoreDual;
	int msg;
	msg=getScorePrimal(&scorePrimal);
	if(msg!=0) { return msg; }
	msg=getScoreDual(&scoreDual);
	if(msg!=0) { return msg; }
	
	(*result)=scorePrimal-scoreDual;

	eprintf("score primal: %e\n",scorePrimal);
	eprintf("score dual: %e\n",scoreDual);
	eprintf("error: %e\n",*result);
	
	if (!std::isfinite((*result))) {
		return MSG_NANINERROR;
	}

	return 0;
	

}


int TSinkhornSolverBarycenterKLMarginals::getScorePrimal(double * const result) {
	(*result)=0;

	// compute all marginals and current barycenter (mean of inner marginals)
	std::vector<TMarginalVector> innerMarginals(nMarginals);
	std::vector<TMarginalVector> outerMarginals(nMarginals);
	TMarginalVector barycenter=TMarginalVector::Constant(zres,0.);
	for(int i=0;i<nMarginals;i++) {
		innerMarginals[i]=v[i].cwiseProduct(kernelT[i]*u[i]);
		outerMarginals[i]=u[i].cwiseProduct(kernel[i]*v[i]);
		barycenter+=weights[i]*innerMarginals[i];
	}

	for(int i=0;i<nMarginals;i++) {
		// F_X, first marginal
		(*result)+=kappa*weights[i]*TSinkhornSolverKLMarginals::KL(outerMarginals[i],
				Eigen::Map<TMarginalVector>(mu[i],res[i]),KLTHRESH);
		// F_Y, second marginal
		(*result)+=kappa*weights[i]*TSinkhornSolverKLMarginals::KL(innerMarginals[i],barycenter,KLTHRESH);
		// KL term wrt. kernel (but without constant kernel mass term)
		
		// first: compute effective dual variables
		TMarginalVector alphaEff=Eigen::Map<TMarginalVector>(alpha[i],res[i])+eps*u[i].array().log().matrix();
		TMarginalVector betaEff=Eigen::Map<TMarginalVector>(beta[i],zres)+eps*v[i].array().log().matrix();
		
		(*result)+=weights[i]*(outerMarginals[i].dot(alphaEff)
				+innerMarginals[i].dot(betaEff)
				-eps*outerMarginals[i].sum());
		
	}
	return 0;



}


int TSinkhornSolverBarycenterKLMarginals::getScoreDual(double * const result) {
	(*result)=0;


	for(int i=0;i<nMarginals;i++) {
		// F_X*, first marginal
		TMarginalVector alphaEff=Eigen::Map<TMarginalVector>(alpha[i],res[i])+eps*u[i].array().log().matrix();
		(*result)-=kappa*weights[i]*TSinkhornSolverKLMarginals::KLDual((-1/kappa)*alphaEff,				
				Eigen::Map<TMarginalVector>(mu[i],res[i]));
		// ignore F_Y* dual contribution. it is an indicator function
		// KL term wrt. kernel (but without constant kernel mass term)
		(*result)-=eps*weights[i]*(u[i].dot(kernel[i]*v[i]));
		
	}

	return 0;



}

