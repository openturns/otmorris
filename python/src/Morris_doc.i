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
We note :math:`\cM:\Rset^p \mapsto \Rset^q` with :math:`\cM(\vect{x})= \vect{y}`.

The Morris method is a screening method, which is known to be very efficient in case of huge number of input parameters (p >> 1).
It is a qualitative sensitivity analysis method which is based on design of experiments and allows to identify the few important factors at a cost of r * (p + 1) simulations.
The experiments are of type OAT (One At Time); i.e. only one parameter vary at a time.

The method helps to split input parameters into three groups:

 - Those with negligible effects on the output,
 - Those with significant and linear effects on the output,
 - Those with significant and non linear (or with interactions) effects on the
   output.

 The method rely on input designs defined in the hypersphere unit. To sum up the key points of the method, we consider a point named :math:`\vect{x^*}` in this hypersphere and a parameter :math:`\delta` (parameter of discretization if we consider a regular experiment for example). Starting from the point, we choose randomly one direction by increasing\slash decreasing one component one component of the point :math:`\vect{x^*}` with :math:`\delta`. Conditionnaly to this direction, we choose then the p-1 directions by randomly selecting one direction at time. We get then a trajectory (path).

The Morris method rely on the evaluation of elementary effects which are defined as follow:

.. math::

    d_{i}(\vect{x}^k) = \frac{\cM(x_1^k,\hdots, x_{i-1}^k, x_i^k + \delta,\hdots, x_p^k) - \cM(x_1^k,\hdots, x_{i-1}^k, x_i^k,\hdots, x_p^k)}{\delta}

With :math:`N` trajectories, we get the mean and standard deviation of these effects (we consider the mean of absolute mean effects in our case). The mean explains the sensitivity wheras the standard deviation explains the interactions and non linear effects.

With the first constructor, we consider that input experiment has been generated thanks to the :class:`~otmorris.MorrisExperiment` and output is evaluated outside the platform.
With second constructor, the output is evaluated inside the platform.

Examples
--------
>>> import openturns as ot
>>> import otmorris
>>> # Define model
>>> ot.RandomGenerator.SetSeed(1)
>>> alpha = ot.DistFunc.rNormal(10)
>>> beta = ot.DistFunc.rNormal(84)
>>> gamma = ot.DistFunc.rNormal(280)
>>> b0 = ot.DistFunc.rNormal()
>>> model = otmorris.MorrisFunction(alpha, beta, gamma, b0)
>>> # Number of trajectories
>>> r = 5
>>> # Define a k-grid level (so delta = 1/(k-1))
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
marginal : int
    Output marginal of interest

Returns
-------
mean: :py:class:`openturns.Point`
    The mean effects.
"

// ---------------------------------------------------------------------

%feature("docstring") OTMORRIS::Morris::getMeanAbsoluteElementaryEffects
"Get the mean of absolute elementary effects.

Parameters
----------
marginal : int
    Output marginal of interest

Returns
-------
mean: :py:class:`openturns.Point`
    The mean effects.
"

// ---------------------------------------------------------------------

%feature("docstring") OTMORRIS::Morris::getStandardDeviationElementaryEffects
"Get the standard deviation of elementary effects.

Parameters
----------
marginal : int
    Output marginal of interest

Returns
-------
mean: :py:class:`openturns.Point`
    The standard effects
"

// ---------------------------------------------------------------------

%feature("docstring") OTMORRIS::Morris::getInputSample
"Accessor to the input sample.

Returns
-------
inputSample : :py:class:`openturns.Sample`
    The input sample
"

// ---------------------------------------------------------------------

%feature("docstring") OTMORRIS::Morris::getOutputSample
"Accessor to the output sample.

Returns
-------
inputSample : :py:class:`openturns.Sample`
    The output sample
"

// ---------------------------------------------------------------------

%feature("docstring") OTMORRIS::Morris::drawElementaryEffects
"Draw elementary effects.

Plots mean vs standard deviation of elementary effects.

Parameters
----------
marginal : int
    Output marginal of interest
absoluteMean : bool, default=True
    Whether to use absolute mean

Returns
-------
mean: :py:class:`openturns.Graph`
    The elementary effects graph
"
