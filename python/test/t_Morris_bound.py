#!/usr/bin/env python

import openturns as ot
import openturns.testing as ott
import otmorris

ot.RandomGenerator.SetSeed(0)

poutre = ot.SymbolicFunction(
    ["L", "b", "h", "E", "F"], ["F * L^3 / (48 * E * b * h^3 / 12)"]
)

L = ot.ParametrizedDistribution(ot.LogNormalMuSigmaOverMu(5.0, 0.02))
b = ot.ParametrizedDistribution(ot.LogNormalMuSigmaOverMu(2.0, 0.05))
h = ot.ParametrizedDistribution(ot.LogNormalMuSigmaOverMu(0.4, 0.05))
E = ot.ParametrizedDistribution(ot.LogNormalMuSigmaOverMu(3e4, 0.12))
F = ot.ParametrizedDistribution(ot.LogNormalMuSigmaOverMu(0.1, 0.20))
distribution = ot.JointDistribution([L, b, h, E, F])
dim = distribution.getDimension()

level_number = 4
trajectories = 10
jump_step = int(level_number / 2)

bounds = ot.Interval(
    [marginal.computeQuantile(0.01)[0] for marginal in [L, b, h, E, F]],
    [marginal.computeQuantile(0.99)[0] for marginal in [L, b, h, E, F]],
)
experiment = otmorris.MorrisExperimentGrid([level_number] * dim, trajectories, bounds)
experiment.setJumpStep(ot.Indices([jump_step] * dim))

input_sample = experiment.generate()
assert input_sample.getSize() == trajectories * (dim + 1)
output_sample = poutre(input_sample)

morris = otmorris.Morris(input_sample, output_sample, bounds)
ott.assert_almost_equal(
    morris.getMeanAbsoluteElementaryEffects(),
    [0.0002595408194356193, 0.0002183769445646401, 0.0006299108455680803,
     0.0005180755623338556, 0.0008016325658652563],
)
ott.assert_almost_equal(
    morris.getMeanElementaryEffects(),
    [0.0002595408194356193, -0.0002183769445646401, -0.0006299108455680803,
     -0.0005180755623338556, 0.0008016325658652563],
)
ott.assert_almost_equal(
    morris.getStandardDeviationElementaryEffects(),
    [0.00010838674923520657, 0.00010581209849382817, 0.0002519423028317143,
     0.00016306903214534856, 0.00028396268061265577],
)
ee = morris.getElementaryEffects()
assert ee.getSize() == trajectories
assert ee.getDimension() == dim
ott.assert_almost_equal(ee.computeMean(), morris.getMeanElementaryEffects())
