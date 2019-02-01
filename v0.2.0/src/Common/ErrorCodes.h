#ifndef ErrorCodes_H_
#define ErrorCodes_H_

// fundamental implementation errors

const int ERR_BASE_NOTIMPLEMENTED=11; // something is not yet implemented. e.g. called empty virtual placeholder method of base class.
const int ERR_BASE_GENERIC=12; // generic error code for something which does not have a unique code.
const int ERR_BASE_UNINIT=13; // something was not properly initialized upon calling
const int ERR_BASE_NOSOL=14; // LP subsolver failed to solve a problem instance

const int ERR_BASE_VARLIST_MERGE_RESFAIL=101; // var list merge failed due to different res values


// measure and marginal preprocessing errors

const int ERR_PREP_INIT_MUXNEG=10001; // minimum value of muX is not strictly positive (short cut solver needs strictly positive masses everywhere)
const int ERR_PREP_INIT_MUYNEG=10002; // same for muY
const int ERR_PREP_INIT_POINTCLOUD2D=10003; // marginal point clouds must be 2d arays

const int ERR_PREP_INIT_DIMMISMATCH=10004; // dimensions of marginals do not coincide. posX and posY are point clouds with in different dimensions

const int ERR_PREP_TRUNC_MUXNEG=10101; // same as ERR_PREP_INIT_MUXNEG, but after mass truncation for integer solvers has been applied,
				      // to separate sources of error					
const int ERR_PREP_TRUNC_MUYNEG=10102; // same for muY

// shortcut solver
const int ERR_SHORTCUT_SUPPORTROWEMPTY=210001; // extracting of support for constructing shielding neighbourhood yielded empty row

// multi scale method
const int ERR_MULTISCALE_SHORTCUTUNSOLVED=20001; // shortcut solver did not succeed.
const int ERR_MULTISCALE_EXCEEDEDLEVELS=20002; // set layer above allowed range
const int ERR_MULTISCALE_SUBSOLVERNOTINITIALIZED=20003; // instance of TShortCutSubSolverInterface is not properly initialized upon calling a method

// basis refinement algorithm

const int ERR_MULTISCALE_BASISREFINEMENT_MARGX=20111; // refined feasible basis coupling has wrong x marginal
const int ERR_MULTISCALE_BASISREFINEMENT_MARGY=20112; // refined feasible basis coupling has wrong y marginal
const int ERR_MULTISCALE_BASISREFINEMENT_BASISENTRIES=20113; // refined feasible basis has wrong number of non-zero entries
const int ERR_MULTISCALE_BASISREFINEMENT_DEPLETED=20114; // mass mismatch when filling refined coupling (too much mass)
const int ERR_MULTISCALE_BASISREFINEMENT_EXCESS=20115; // mass mismatch when filling refined coupling (too little mass)

const int ERR_MULTISCALE_BASISREFINEMENT_COARSEBASISENTRIES=20116; // coarse input basis has wrong number of non-zero entries

// CPLEX sub solver
const int ERR_CPLEX_BASISEXPORT_BASISENTRIES=14001; // number of entries in extracted basis info was wrong

// Lemon sub solver
const int ERR_LEMON_RSLT_ERROR_NOSOL=15001;
const int ERR_LEMON_RSLT_ERROR_UNBOUNDED=15002;
const int ERR_LEMON_RSLT_ERROR_UNDEFINED=15003;

// lp_solve sub solver
const int ERR_LPSOLVE_BASISEXPORT_SIGN=16001; // sign of basis variable data was positive (should all be negative)
const int ERR_LPSOLVE_BASISEXPORT_COUNT=16002; // nr of basis variables not right (should be xres+yres-1)
const int ERR_LPSOLVE_BASISIMPORT_COUNT=16003; // nr of basis variables not right (should be xres+yres-1)

// sparse simplex solver

const int ERR_SPARSESIMPLEX_TRUNCATION=13001; // formal truncation (double -> int casting) before calling subsolver caused different masses

// epsilon scaling
const int ERR_EPSSCALING_EMPTYCOARSESUBLIST=50001; // empty eps sub list on coarsest hierarchy level

#endif
