// SWIG file

%{
#include "otmorris/Morris.hxx"
%}

%include Morris_doc.i

%copyctor OTMORRIS::Morris;

%include otmorris/Morris.hxx

%pythoncode %{

import openturns as ot
import numpy as np

class MorrisFunction(ot.OpenTURNSPythonFunction):
    """
    Morris test function for sensitivity analysis.

    This function has input dimension 20 and output dimension 1.

    References
    ----------
    - M. D. Morris, 1991, Factorial sampling plans for preliminary
      computational experiments, Technometrics, 33, 161-174.

    Examples
    --------
    >>> import openturns as ot
    >>> ot.RandomGenerator.SetSeed(123)
    >>> b0_random = ot.DistFunc.rNormal()
    >>> b1_random = ot.DistFunc.rNormal(10)
    >>> b2_random = ot.DistFunc.rNormal(175)
    >>> morrisFunction = ot.Function(MorrisFunction(b0_random, b1_random, b2_random))
    >>> dimension = morrisFunction.getInputDimension()
    >>> distribution = ot.ComposedDistribution([ot.Uniform(0.0, 1.0)] * dimension)
    >>> input_sample = distribution.getSample(10)
    >>> output_sample = morrisFunction(input_sample)
    """

    def __init__(self, b0_random=0.0, b1_random=ot.Point(10), b2_random=ot.Point(175)):
        """
        Create the Morris function.

        Parameters
        ----------
        b0_random : float, optional
            The constant term. Default is 0.0.
        b1_random : ot.Point(10), optional
            Random linear coefficients for dimensions 11-20. Default is zeros.
        b2_random : ot.Point(175), optional
            Random quadratic coefficients. Default is zeros.
        """
        super().__init__(20, 1)
        self.b0_random = float(b0_random)

        # Initialize b1_random (default to zeros)
        b1_random = np.array(b1_random)
        assert len(b1_random) == 10, f"b1_random must have length 10, got {len(b1_random)}"
        self.b1 = np.concatenate([np.full(10, 20.0), b1_random])

        # Initialize b2_random (default to zeros)
        b2_random = np.array(b2_random)
        assert len(b2_random) == 175, f"b2_random must have length 175, got {len(b2_random)}."

        # Precompute b2 matrix (20x20) using numpy for speed
        self.b2 = np.zeros((20, 20))

        # Fill quadratic terms (i < j)
        random_index = 0
        for i in range(20):
            for j in range(i + 1, 20):
                if i < 6 and j < 6:
                    # The first 15 terms are deterministic
                    self.b2[i, j] = -15.0
                else:
                    # The other 175 coefficients are random
                    self.b2[i, j] = b2_random[random_index]
                    random_index += 1

        # Precompute quadratic indices
        self.quad_i, self.quad_j = np.triu_indices(20, k=1)
        self.b2_upper = self.b2[self.quad_i, self.quad_j]

        # Precompute b3 tensor (20x20x20) - sparse, only 5x5x5 block is non-zero
        self.b3 = np.zeros((20, 20, 20))
        self.b3[:5, :5, :5] = -10.0

        # Precompute cubic term indices (i < j < k)
        cubic_indices = []
        cubic_coeffs = []
        for i in range(20):
            for j in range(i + 1, 20):
                for k in range(j + 1, 20):
                    coeff = self.b3[i, j, k]
                    if coeff != 0:
                        cubic_indices.append((i, j, k))
                        cubic_coeffs.append(coeff)
        self.cubic_indices = np.array(cubic_indices)
        self.cubic_coeffs = np.array(cubic_coeffs)

        # Precompute b4 tensor (20x20x20x20) - sparse, only 4x4x4x4 block is non-zero
        # Store only non-zero quartic terms
        quartic_indices = []
        quartic_coeffs = []
        for i in range(4):
            for j in range(i + 1, 4):
                for k in range(j + 1, 4):
                    for ell in range(k + 1, 4):
                        quartic_indices.append((i, j, k, ell))
                        quartic_coeffs.append(5.0)
        self.quartic_indices = np.array(quartic_indices)
        self.quartic_coeffs = np.array(quartic_coeffs)

    def _exec(self, x):
        """Evaluate the Morris function at point x (vectorized)."""
        # Convert to numpy array
        X = np.array(x)

        # Transform to w coordinates (vectorized)
        w = (X - 0.5) * 2.0

        # Apply nonlinear transformation to indices 2, 4, 6
        for k in [2, 4, 6]:
            w[k] = 2.0 * (1.1 * X[k] / (X[k] + 0.1) - 0.5)

        # Start with constant term
        y = self.b0_random

        # Add linear terms (vectorized dot product)
        y += np.dot(w, self.b1)

        # Add quadratic terms (vectorized using precomputed indices)
        y += np.sum(self.b2_upper * w[self.quad_i] * w[self.quad_j])

        # Add cubic terms (vectorized using precomputed non-zero terms)
        w_cubic = (
            w[self.cubic_indices[:, 0]]
            * w[self.cubic_indices[:, 1]]
            * w[self.cubic_indices[:, 2]]
        )
        y += np.sum(self.cubic_coeffs * w_cubic)

        # Add quartic terms (vectorized using precomputed non-zero terms)
        w_quartic = (
            w[self.quartic_indices[:, 0]]
            * w[self.quartic_indices[:, 1]]
            * w[self.quartic_indices[:, 2]]
            * w[self.quartic_indices[:, 3]]
        )
        y += np.sum(self.quartic_coeffs * w_quartic)

        return [float(y)]

%}
