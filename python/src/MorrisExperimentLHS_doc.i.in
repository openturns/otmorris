%feature("docstring") OTMORRIS::MorrisExperimentLHS
"MorrisExperimentLHS builds experiments for the Morris method using a centered LHS design as input starting.

Parameters
----------
lhsDesign : :py:class:`openturns.Sample`
    Initial design
N : int
    Number of trajectories
bounds : :py:class:`openturns.Interval`, optional
    Bounds of the domain, by default it is defined on :math:`[0,1]^d`.

Notes
-----
The lhs sample must be centered, ie from :py:class:`openturns.LHSExperiment` with randomShift=False.

The method consists in generating trajectories (paths) by randomly selecting their initial points from the lhs design.
If number of trajectories is lesser than the lhsDesign's size, we enforce the selection of the starting point using
:py:class:`openturns.KPermutationsDistribution` which ensure full different trajectories.

Examples
--------
>>> import openturns as ot
>>> import otmorris
>>> ot.RandomGenerator.SetSeed(1)
>>> r = 5
>>> # Define experiments in [0,1]^2
>>> size = 20
>>> # Generate an LHS design
>>> dist = ot.ComposedDistribution([ot.Uniform(0, 1)] * 2)
>>> # should be centered so randomShift=False
>>> lhs_experiment = ot.LHSExperiment(dist, size, True, False)
>>> lhsDesign = lhs_experiment.generate()
>>> experiment = otmorris.MorrisExperimentLHS(lhsDesign, r)
>>> X = experiment.generate()
"
