%feature("docstring") OTMORRIS::MorrisExperiment
"Base class for the Morris method experiments.

See also
--------
MorrisExperimentGrid, MorrisExperimentLHS"

// ---------------------------------------------------------------------

%feature("docstring") OTMORRIS::MorrisExperiment::getBounds
"Get the bounds of the domain.

Returns
-------
bounds : :py:class:`openturns.Interval`
    Bounds of the domain, default is :math:`[0,1]^d`
"

// ---------------------------------------------------------------------

%feature("docstring") OTMORRIS::MorrisExperiment::generate
R"RAW(Generate points according to the type of the experiment.

Returns
-------
sample : :py:class:`openturns.Sample`
    Points that constitute the design of experiment, of size :math:`N \times (p+1)`
)RAW"
