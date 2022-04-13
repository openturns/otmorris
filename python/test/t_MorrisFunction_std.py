#!/usr/bin/env python

import openturns as ot
import openturns.testing as ott
import otmorris
import time

f = ot.Function(otmorris.MorrisFunction())
dim = f.getInputDimension()
x1 = [0.3] * dim
x2 = [0.4] * dim
x3 = [0.5] * dim
y1 = f(x1)
y2 = f(x2)
y3 = f(x3)
print(y1, y2, y3)
ott.assert_almost_equal(y1, [-17.738])
ott.assert_almost_equal(y2, [19.2912])
ott.assert_almost_equal(y3, [50])

X = ot.Normal([0.4] * dim, [0.1] * dim, ot.CorrelationMatrix(dim))
N = 1000
x = X.getSample(N)
t0 = time.time()
y = f(x)
t1 = time.time()
print(N / (t1 - t0), 'eval/s')
