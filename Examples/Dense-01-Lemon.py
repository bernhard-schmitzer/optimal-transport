# load libraries
from lib.header_script import *

# assume that CPLEX back-end is installed
# import dense CPLEX solver interface
import Solvers.OT_Lemon as OT_Lemon

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

# the LEMON solver internally only works with integer marginals and cost function entries.
# therefore marginals and cost must be truncated to a discrete grid, for which we must choose resolutions.
# this truncation usually inflicts a small error in the original problem, which can however be made negligibly small

# truncating measures to dense units. truncationMeasures gives the size of each mass atom.
#	this means that the total mass of both marginals will actually be of the order of 1/truncationMeasures,
#	which has to be taken into account when looking at the final primal and dual score

truncationMeasures=1E-9
OTTools.TruncateMeasures(muX,muY,truncationMeasures)

# the variable truncationCost determines the resolution of the truncation grid for the cost function
truncationCost=1E-3


# compute dense cost function
p=2. # exponent for Euclidean distance
c=OTTools.getEuclideanCostFunction(posX,posY,p=p)




print("Solving dense.",flush=True)
time1=datetime.datetime.now()
# call dense solver from OT_Lemon module. cScale is the parameter that handles the truncation of the cost function
# algorithm=0 uses the Lemon Network Simplex algorithm
# algorithm=1 uses the Cost Scaling algorithm
result=OT_Lemon.solveDense(c,muX,muY,cScale=truncationCost,algorithm=0)
time2=datetime.datetime.now()
deltaTDense=(time2-time1).total_seconds()
print("\trequired time: ",time2-time1)

# OT_Lemon.solveDense returns a dictionary with keys [mu,alpha,beta,status], containing the optimal coupling,
#	the optimal dual variables and the status of the solver (which is 0 if everything went fine).
# Just for illustration use primal and dual optimizers to compute primal and dual optimal score.
# Note that the scores will be basically scaled by 1/truncationMeasures
print("\tprimal score: ", np.sum(result["mu"]*c))
print("\tdual score: ", np.sum(result["alpha"]*muX)+np.sum(result["beta"]*muY) )
