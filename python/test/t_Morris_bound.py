#!/usr/bin/env python

import openturns as ot
import otmorris

poutre = ot.SymbolicFunction(
    ["L", "b", "h", "E", "F"], ["F * L^3 / (48 * E * b * h^3 / 12)"]
)

# define the model
L = ot.ParametrizedDistribution(ot.LogNormalMuSigmaOverMu(5.0, 0.02))
b = ot.ParametrizedDistribution(ot.LogNormalMuSigmaOverMu(2.0, 0.05))
h = ot.ParametrizedDistribution(ot.LogNormalMuSigmaOverMu(0.4, 0.05))
E = ot.ParametrizedDistribution(ot.LogNormalMuSigmaOverMu(3e4, 0.12))
F = ot.ParametrizedDistribution(ot.LogNormalMuSigmaOverMu(0.1, 0.20))
list_marginals = [L, b, h, E, F]
distribution = ot.ComposedDistribution(list_marginals)
distribution.setDescription(("L", "b", "h", "E", "F"))
dim = distribution.getDimension()

level_number = 4
trajectories = 10
jump_step = int(level_number / 2)
levels = [level_number] * dim

# set the bounds of the grid experiment
bound = ot.Interval(
    [marginal.computeQuantile(0.01)[0] for marginal in list_marginals],
    [marginal.computeQuantile(0.99)[0] for marginal in list_marginals],
)
experiment = otmorris.MorrisExperimentGrid(levels, bound, trajectories)
experiment.setJumpStep(ot.Indices([jump_step] * dim))

# create and compute the design of experiments
input_sample = experiment.generate()
print("Morris experiment generated from grid = ")
print(input_sample)
output_sample = poutre(input_sample)

# run the Morris analysis
morris = otmorris.Morris(input_sample, output_sample, bound)
print("E(|EE|)  = ", morris.getMeanAbsoluteElementaryEffects())
print("E(EE)  = ", morris.getMeanElementaryEffects())
print("V(|EE|)^{1/2} = ", morris.getStandardDeviationElementaryEffects())
