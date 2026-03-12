#!/usr/bin/env python

import openturns as ot
import openturns.testing as ott
import otmorris
import time
import numpy as np

ot.RandomGenerator.SetSeed(1)
b0_random = ot.DistFunc.rNormal()
b1_random = ot.DistFunc.rNormal(10)
b2_random = ot.DistFunc.rNormal(175)
g = ot.Function(otmorris.MorrisFunction(b0_random, b1_random, b2_random))
dim = g.getInputDimension()

# Check accuracy
x = np.linspace(0.0, 1.0, 20)
y = g(x)
ott.assert_almost_equal(y, [-65.75761172072895])

# Check speed
X = ot.ComposedDistribution([ot.Uniform(0.0, 1.0)] * 20)
N = 1000
input_sample = X.getSample(N)
t0 = time.time()
output_sample = g(input_sample)
t1 = time.time()
elapsed_time = t1 - t0
print(f"{N / elapsed_time} eval/s")
