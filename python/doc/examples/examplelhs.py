import openturns as ot
import otmorris
from otmorris.plot_sensitivity import PlotEE

# Number of trajectories
r = 10
dim = 20
# Define an LHS experiment of size 50 in [0, 1]^20
size = 50
dist = ot.ComposedDistribution([ot.Uniform(0, 1)] * dim)
lhs_experiment = ot.LHSExperiment(dist, size, True, False)
lhsDesign = lhs_experiment.generate()
morris_experiment = otmorris.MorrisExperimentLHS(lhsDesign, r)
bounds = ot.Interval(dim) # [0, 1]^2
X = morris_experiment.generate()
f = ot.Function(otmorris.MorrisFunction())
Y = f(X)
# Evaluate Elementary effects (ee)
ee = otmorris.Morris(X, Y, bounds)
# Compute mu/sigma
mean = ee.getMeanAbsoluteElementaryEffects()
sigma = ee.getStandardDeviationElementaryEffects()
fig = PlotEE(ee, title="Elementary Effects using LHS")
fig.show()
