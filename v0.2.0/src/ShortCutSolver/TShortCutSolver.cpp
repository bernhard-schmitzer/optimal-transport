#include"TShortCutSolver.h"


TShortCutSolver::TShortCutSolver(
		TShortCutCouplingHandlerInterface *_couplingHandler, TShortCutSubSolverInterfaceBase *_solverInterface,
		TShieldGeneratorBase *_shieldGenerator, int _violationCheckMethod) {
	solverInterface=_solverInterface;
	couplingHandler=_couplingHandler;
	shieldGenerator=_shieldGenerator;
	violationCheckMethod=_violationCheckMethod;
	report.objective=0;
	report.solved=0;
	report.steps=0;

	report.t_solving=0;
	report.t_shielding=0;

}

TShortCutSolver::~TShortCutSolver() {
}

int TShortCutSolver::initialize() {
	report.objective=0;
	report.solved=0;
	report.steps=0;

	report.t_solving=0;
	report.t_shielding=0;

	report.n_variables=0;
	report.n_shielding_queries=0;

	return 0;
}

int TShortCutSolver::solve() {
	int msg;
	msg=solverInterface->solve();
	return msg;
}

double TShortCutSolver::getObjective() {
	return solverInterface->getObjective();
}

int TShortCutSolver::step() {
	return step(1);
}

int TShortCutSolver::step(int steps) {
	int n;
	bool violation;
	int msg;
	double oldObjective;
	TVarListHandler *newXVars;
	int time1,time2;

	n=0;
	violation=true;
	msg=0;
	while((n<steps)&&(violation)&&(msg==0)) {
		
		oldObjective=getObjective();
		time1=clock();
		msg=solverInterface->solve();
		time2=clock();
		report.t_solving+=(time2-time1);
		if(msg!=0) {
			// something low level went wrong during solving
			return msg;
		}
		report.steps++;
		report.objective=getObjective();
		

		// primal convergence check can be done before building shield
		switch(violationCheckMethod) {
		case VCHECK_PRIMAL:
			if((getObjective()>=oldObjective) && (n>0)) {
				violation=false;
			}
			break;
		}

		if(violation) {
			newXVars=couplingHandler->getSupport();
			
			// check that each row has at least one entry
			if(VarListTools_HasEmptyRows(newXVars)) {
				return ERR_SHORTCUT_SUPPORTROWEMPTY;
			}

			time1=clock();
			shieldGenerator->generateShield(newXVars,newXVars);
			time2=clock();
			report.t_shielding+=(time2-time1);
			report.n_variables=newXVars->total;
			report.n_shielding_queries=shieldGenerator->n_shielding_queries;

			newXVars->sort();

			// subsolver may extract status info from current state and add some variables to newXVars if desired
			solverInterface->prepareUpdate(newXVars);

			couplingHandler->updateXVars(newXVars,false);


			// dual convergence checks come after constructing the shield
			switch(violationCheckMethod) {
			case VCHECK_DUAL:
				violation=couplingHandler->dualViolationCheck(false, solverInterface->alpha, solverInterface->beta);
				break;
			case VCHECK_DUAL_PROJECT:
				violation=couplingHandler->dualViolationCheck(true, solverInterface->alpha, solverInterface->beta);
				break;
			}
			n++;
			
			// subsolver updates to new xVars setting
			msg=solverInterface->update();
			if(msg!=0) {
				// something went wrong during updating solver
				return msg;
			}

		}
	}


	if(msg!=0) {
		return msg;
	}

	if(!violation) {
		report.solved=1;
		return 0;
	} else {
		report.solved=0;
		return 1;
	}

}


