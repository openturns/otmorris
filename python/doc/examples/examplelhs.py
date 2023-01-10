import openturns as ot
import otmorris
from openturns.viewer import View

# use the reference 20-d function from the Morris paper
f = ot.Function(otmorris.MorrisFunction())
dim = f.getInputDimension()

# Number of trajectories
r = 10

# Define an LHS experiment of size 50 in [0, 1]^20
size = 50
dist = ot.ComposedDistribution([ot.Uniform(0, 1)] * dim)
lhs_experiment = ot.LHSExperiment(dist, size, True, False)
lhsDesign = lhs_experiment.generate()
morris_experiment = otmorris.MorrisExperimentLHS(lhsDesign, r)
bounds = ot.Interval(dim) # [0, 1]^2
X = morris_experiment.generate()
Y = f(X)

# Evaluate Elementary effects (ee)
morris = otmorris.Morris(X, Y, bounds)

# Compute mu/sigma
mean = morris.getMeanAbsoluteElementaryEffects()
sigma = morris.getStandardDeviationElementaryEffects()
graph = morris.drawElementaryEffects(0)
graph.setTitle("Elementary Effects using LHS")
View(graph).show()
