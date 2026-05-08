#!/usr/bin/env python

import openturns as ot
import openturns.testing as ott
import otmorris
import time


def linspace(xmin, xmax, npoints):
    """Returns a sample created from a regular grid
    from xmin to xmax with npoints points."""
    step = (xmax - xmin) / (npoints - 1)
    rg = ot.RegularGrid(xmin, step, npoints)
    vertices = rg.getVertices()
    return vertices.asPoint()


ot.RandomGenerator.SetSeed(1)
b0Random = ot.DistFunc.rNormal()
b1Random = ot.DistFunc.rNormal(10)
b2Random = ot.DistFunc.rNormal(175)
g = otmorris.BuildMorrisFunction(b0Random, b1Random, b2Random)
dim = g.getInputDimension()

# Check accuracy
x = linspace(0.0, 1.0, 20)
y = g(x)
ott.assert_almost_equal(y, [-65.75761172072895])

# Check speed
X = ot.JointDistribution([ot.Uniform(0.0, 1.0)] * 20)
N = 10000
input_sample = X.getSample(N)
t0 = time.time()
output_sample = g(input_sample)
t1 = time.time()
elapsed_time = t1 - t0
print(f"{N / elapsed_time} eval/s")
