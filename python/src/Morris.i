// SWIG file

%{
#include "otmorris/Morris.hxx"
%}

%include Morris_doc.i

%include otmorris/Morris.hxx
namespace OTMORRIS { %extend Morris { Morris(const Morris & other) { return new OTMORRIS::Morris(other); } } }

%pythoncode %{

import openturns as ot

class MorrisFunction(ot.OpenTURNSPythonFunction):
    """
    The non-monotonic function of Morris f: R^20 -> R

    Reference:
      M. D. Morris, 1991, Factorial sampling plans for preliminary
      computational experiments,Technometrics, 33, 161--174.

    Examples
    --------
    >>> import openturns as ot
    >>> ot.RandomGenerator.SetSeed(123)
    >>> b0 = ot.DistFunc.rNormal()
    >>> alpha = ot.DistFunc.rNormal(10)
    >>> beta =  ot.DistFunc.rNormal(6*14)
    >>> gamma =  ot.DistFunc.rNormal(20*14)
    >>> f = ot.Function(MorrisFunction(alpha, beta, gamma, b0))
    >>> input_sample = ot.ComposedDistribution([ot.Uniform(0,1)] * 20).getSample(20)
    >>> output_sample = f(input_sample)

    """
    def __init__(self, alpha = ot.Point(10), beta = ot.Point(14*6),\
                 gamma = ot.Point(20*14), b0 = 0.0):
        ot.OpenTURNSPythonFunction.__init__(self, 20, 1)
        self.b0 = float(b0)
        # Check alpha dimension
        assert(len(alpha) == 10)
        self.b1 = [20] * 10 + list(alpha)
        # Check beta and gamma dimension
        assert(len(beta) == 6 * 14)
        assert(len(gamma) == 20 * 14)
        self.b2 = [[0] *20] * 20
        for i in range(6):
            for j in range(6):
                self.b2[i][j] = -15.0
        # Take into account beta
        k = 0
        for i in range(6):
            for j in range(14):
                self.b2[i][j + 6] = beta[k]
                k = k + 1
        # Take into account gamma
        k = 0
        for i in range(6, 20):
            for j in range(20):
                self.b2[i][j] = gamma[k]

        # b3
        self.b3 = [[[0]*20]*20]*20
        for i in range(5):
            for j in range(5):
                for k in range(5):
                    self.b3[i][j][k] = -10.0
        # b4
        self.b4 = [[[[0]*20]*20]*20]*20
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for l in range(4):
                        self.b4[i][j][k][l] = 5

    def _exec(self, x):
        assert (len(x)==20)
        b1 = self.b1
        b2 = self.b2
        b3 = self.b3
        b4 = self.b4
        # X is a list, transform it into array
        X = ot.Point(x)
        w = (X - [0.5]*20)*2
        for k in [2,4,6]:
           w[k] = 2.0 * (1.1 * X[k] / (X[k] + 0.1) - 0.5)
        y = w.dot(b1)
        # Morris function
        for i in range(19):
            for j in range(i + 1, 20):
                y +=  b2[i][j] * w[i] * w[j]
        for i in range(18):
            for j in range(i + 1, 19):
                for k in range(j + 1, 20):
                    y += b3[i][j][k] * w[i] * w[j] * w[k]

        for i in range(17):
            for j in range(i + 1, 18):
                for k in range(j + 1, 20):
                    for l in range(k + 1, 20):
                        y += b4[i][j][k][l] * w[i] * w[j] * w[k] * w[l]

        return [y + self.b0]
%}
