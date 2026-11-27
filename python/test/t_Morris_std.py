#!/usr/bin/env python

import openturns as ot
import openturns.testing as ott
import otmorris

ot.RandomGenerator.SetSeed(0)
r = 5
morris_experiment = otmorris.MorrisExperimentGrid([5] * 2, r)
grid_bound = morris_experiment.getBounds()
sample1 = morris_experiment.generate()
assert sample1.getSize() == r * 3
assert sample1.getDimension() == 2

size = 20
dist = ot.JointDistribution(2 * [ot.Uniform(0, 1)])
lhsDesign = ot.LHSExperiment(dist, size, True, False).generate()
morris_experiment_lhs = otmorris.MorrisExperimentLHS(lhsDesign, r)
lhs_bound = morris_experiment_lhs.getBounds()
sample2 = morris_experiment_lhs.generate()
assert sample2.getSize() == r * 3
assert sample2.getDimension() == 2

model = ot.SymbolicFunction(["x", "y"], ["cos(x)*y + sin(y)*x + x*y -0.1"])

morrisEE1 = otmorris.Morris(sample1, model(sample1), grid_bound)
ott.assert_almost_equal(
    morrisEE1.getMeanAbsoluteElementaryEffects(), [0.4495333571083234, 1.4993301542775208]
)
ott.assert_almost_equal(
    morrisEE1.getMeanElementaryEffects(), [0.4495333571083234, 1.4993301542775208]
)
ott.assert_almost_equal(
    morrisEE1.getStandardDeviationElementaryEffects(),
    [0.2817970916125558, 0.320673153854582],
)
ee1 = morrisEE1.getElementaryEffects()
assert ee1.getSize() == r
assert ee1.getDimension() == 2
ott.assert_almost_equal(ee1.computeMean(), morrisEE1.getMeanElementaryEffects())
ott.assert_almost_equal(
    ee1.computeStandardDeviation(), morrisEE1.getStandardDeviationElementaryEffects()
)

morrisEE2 = otmorris.Morris(sample2, model(sample2), lhs_bound)
ott.assert_almost_equal(
    morrisEE2.getMeanAbsoluteElementaryEffects(), [0.7948394944942017, 1.9055894592131428]
)
ott.assert_almost_equal(
    morrisEE2.getStandardDeviationElementaryEffects(),
    [0.5963247385341721, 0.45962545940559757],
)
ee2 = morrisEE2.getElementaryEffects()
assert ee2.getSize() == r
assert ee2.getDimension() == 2
ott.assert_almost_equal(ee2.computeMean(), morrisEE2.getMeanElementaryEffects())
ott.assert_almost_equal(
    ee2.computeStandardDeviation(), morrisEE2.getStandardDeviationElementaryEffects()
)

# overflow check
levels = ot.Indices(168)
levels.fill(5, 0)
morris_experiment = otmorris.MorrisExperimentGrid(levels, r)
