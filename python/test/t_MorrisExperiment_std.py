#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import otmorris

ot.RandomGenerator.SetSeed(0)
r = 5
morris_experiment = otmorris.MorrisExperimentGrid([5, 5], r)
X = morris_experiment.generate()
assert X.getSize() == r * 3
assert X.getDimension() == 2

# Each trajectory step moves a single coordinate by delta=0.25
delta = 0.25
for n in range(r):
    first = n * (X.getDimension() + 1)
    for i in range(X.getDimension()):
        step = X[first + i + 1] - X[first + i]
        nonzero = [abs(v) > 1e-12 for v in step]
        assert sum(nonzero) == 1
        for v in step:
            if abs(v) > 1e-12:
                ott.assert_almost_equal(abs(v), delta)

ot.RandomGenerator.SetSeed(1)
lhs = ot.Sample(
    [
        [0.875, 0.575],
        [0.775, 0.775],
        [0.625, 0.925],
        [0.375, 0.625],
        [0.275, 0.975],
        [0.075, 0.825],
        [0.225, 0.325],
        [0.425, 0.225],
        [0.925, 0.475],
        [0.325, 0.725],
        [0.025, 0.275],
        [0.525, 0.175],
        [0.125, 0.375],
        [0.475, 0.025],
        [0.825, 0.675],
        [0.175, 0.525],
        [0.975, 0.875],
        [0.675, 0.075],
        [0.725, 0.125],
        [0.575, 0.425],
    ]
)

morris_experiment = otmorris.MorrisExperimentLHS(lhs, r)
Y = morris_experiment.generate()
assert Y.getSize() == r * 3
assert Y.getDimension() == 2

# delta = 1 / lhs size = 0.05, single-coordinate steps
for n in range(r):
    first = n * (Y.getDimension() + 1)
    for i in range(Y.getDimension()):
        step = Y[first + i + 1] - Y[first + i]
        nonzero = [abs(v) > 1e-12 for v in step]
        assert sum(nonzero) == 1
        for v in step:
            if abs(v) > 1e-12:
                ott.assert_almost_equal(abs(v), 0.05)
