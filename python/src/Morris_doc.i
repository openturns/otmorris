%feature("docstring") OTMORRIS::Morris
R"RAW(Morris method.

Available constructors:

    Morris(*inputSample, outputSample, interval*)

    Morris(*experiment, model*)

Parameters
----------
inputSample : :py:class:`openturns.Sample`
    Experiment generated thanks to the `generate` method of the :class:`~otmorris.MorrisExperiment`
outputSample : :py:class:`openturns.Sample`
    Response model applied on `inputSample`
interval : :py:class:`openturns.Interval`
    Bounds of the experiment inputs.
experiment : :py:class:`otmorris.MorrisExperiment`
    Morris experiment
model : :py:class:`openturns.Function`
    Response model to be applied on input data

Notes
-----
We note :math:`\model:\Rset^{\inputDim} \mapsto \Rset^{\outputDim}` the physical
model with :math:`\model(\vect{x}) = \vect{y}`.

The Morris method is a screening method, which is known to be very efficient in
case of large number of input parameters (:math:`\inputDim \gg 1`).
It is a qualitative sensitivity analysis method which is based on design of
experiments and allows to identify the few important factors at a cost of
:math:`r * (\inputDim + 1)` simulations.
The experiments are of type OAT (One At Time); i.e. only one parameter varies
at a time.

The method helps to split input parameters into three groups:

 - Those with negligible effects on the output,
 - Those with significant and linear effects on the output,
 - Those with significant and non linear (or with interactions) effects on the
   output.

The method relies on input designs defined in the hypersphere unit.
To sum up the key points of the method, we consider a point named
:math:`\vect{x^*}` in this hypersphere and a parameter :math:`\delta`
(parameter of discretization if we consider a regular experiment for example).
Starting from the point, we choose randomly one direction by increasing or
decreasing one component of the point :math:`\vect{x^*}` with :math:`\delta`.
Conditionally to this direction, we choose then the :math:`\inputDim - 1`
directions by randomly selecting one direction at a time.
We then get a trajectory (path).

The Morris method is based on the evaluation of elementary effects which are
defined as follows:

.. math::

    d_{i}(\vect{x}^k) 
    = \frac{\model(x_1^k, \ldots, x_{i - 1}^k, x_i^k + \delta, \ldots, x_{\inputDim}^k) 
        - \model(x_1^k, \ldots, x_{i - 1}^k, x_i^k, \ldots, x_{\inputDim}^k)}{\delta}.

With N trajectories, we get the mean and standard deviation of these effects
(we consider the mean of absolute mean effects in our case).
The mean explains the sensitivity whereas the standard deviation explains the
interactions and non linear effects.

With the first constructor, we consider that input experiment has been generated
thanks to the :class:`~otmorris.MorrisExperiment` and output is evaluated
outside the platform.
With the second constructor, the output is evaluated inside the platform.

Examples
--------
>>> import openturns as ot
>>> import otmorris
>>> # Define model
>>> ot.RandomGenerator.SetSeed(1)
>>> b0_random = ot.DistFunc.rNormal()
>>> b1_random = ot.DistFunc.rNormal(10)
>>> b2_random = ot.DistFunc.rNormal(175)
>>> model = otmorris.MorrisFunction(b0_random, b1_random, b2_random)
>>> # Number of trajectories
>>> r = 5
>>> # Define a k-grid level (so delta = 1 / (k - 1))
>>> k = 5
>>> morris_experiment = otmorris.MorrisExperimentGrid([k] * 20, r)
>>> X = morris_experiment.generate()
>>> # Evaluation of the model on the design: evaluation outside OT
>>> Y = model(X)
>>> # need the bounds when using X and Y
>>> bounds = morris_experiment.getBounds()
>>> # Evaluation of Morris effects
>>> morris = otmorris.Morris(X, Y, bounds)
>>> # Get mean/sigma effects
>>> mean_effects = morris.getMeanElementaryEffects()
>>> mean_abs_effects = morris.getMeanAbsoluteElementaryEffects()
>>> sigma_effects = morris.getStandardDeviationElementaryEffects()
)RAW"

// ---------------------------------------------------------------------

%feature("docstring") OTMORRIS::Morris::getMeanElementaryEffects
"Get the mean of elementary effects.

Parameters
----------
marginal : int, optional
    Output marginal of interest. Default is 0.

Returns
-------
mean : :py:class:`openturns.Point`
    The mean of elementary effects.
"

// ---------------------------------------------------------------------

%feature("docstring") OTMORRIS::Morris::getMeanAbsoluteElementaryEffects
"Get the mean of absolute elementary effects.

Parameters
----------
marginal : int, optional
    Output marginal of interest. Default is 0.

Returns
-------
meanOfAbsolute : :py:class:`openturns.Point`
    The mean of absolute elementary effects.
"

// ---------------------------------------------------------------------

%feature("docstring") OTMORRIS::Morris::getStandardDeviationElementaryEffects
"Get the standard deviation of the elementary effects.

Parameters
----------
marginal : int, optional
    Output marginal of interest. Default is 0.

Returns
-------
standardDeviation : :py:class:`openturns.Point`
    The standard deviation of the elementary effects.
"

// ---------------------------------------------------------------------

%feature("docstring") OTMORRIS::Morris::getInputSample
"Accessor to the input sample.

Returns
-------
inputSample : :py:class:`openturns.Sample`
    The input sample.
"

// ---------------------------------------------------------------------

%feature("docstring") OTMORRIS::Morris::getOutputSample
"Accessor to the output sample.

Returns
-------
outputSample : :py:class:`openturns.Sample`
    The output sample.
"

// ---------------------------------------------------------------------

%feature("docstring") OTMORRIS::Morris::drawElementaryEffects
"Draw elementary effects.

Plots mean vs standard deviation of elementary effects.

Parameters
----------
marginal : int, optional
    Output marginal of interest. Default is 0.
absoluteMean : bool, optional
    Whether to use absolute mean. Default is True.

Returns
-------
graph : :py:class:`openturns.Graph`
    The elementary effects graph.
"
