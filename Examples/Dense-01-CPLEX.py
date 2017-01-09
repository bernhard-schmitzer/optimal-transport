# load libraries
from lib.header_script import *

# assume that CPLEX back-end is installed
# import dense CPLEX solver interface
import Solvers.OT_CPLEX as OT_CPLEX

########################################
# Problem Setup
########################################
print("Problem setup.",flush=True)

# load two 64x64 images
imgX=sciio.loadmat("data/density-img/f-000-064.mat")["a"]
imgY=sciio.loadmat("data/density-img/f-001-064.mat")["a"]

# preprocessing (add small constant background mass, normalize masses, extract geometric positions of pixels)
(muX,posX)=OTTools.processDensity_Grid(imgX,totalMass=1.,constOffset=1E-6)
(muY,posY)=OTTools.processDensity_Grid(imgY,totalMass=1.,constOffset=1E-6)


# compute dense cost function
p=2. # exponent for Euclidean distance
c=OTTools.getEuclideanCostFunction(posX,posY,p=p)


print("Solving dense.",flush=True)
time1=datetime.datetime.now()
# call dense solver from OT_CPLEX module. getDualVariables=True also retrieves the optimal dual variables
resultDense_CPLEX=OT_CPLEX.solveDense(c,muX,muY,getDualVariables=True)
time2=datetime.datetime.now()
deltaTDense=(time2-time1).total_seconds()
print("\trequired time: ",time2-time1)

# OT_CPLEX.solveDense returns a dictionary with keys [mu,alpha,beta,status], containing the optimal coupling,
#	the optimal dual variables and the status of the solver (which is 0 if everything went fine).
#	alpha and beta are only returned, if OT_CPLEX.solveDense is called with getDualVariables=True.
# Just for illustration use primal and dual optimizers to compute primal and dual optimal score.
print("\tprimal score: ", np.sum(resultDense_CPLEX["mu"]*c))
print("\tdual score: ", np.sum(resultDense_CPLEX["alpha"]*muX)+np.sum(resultDense_CPLEX["beta"]*muY) )
