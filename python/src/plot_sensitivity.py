"""
Plot Morris elementary effects
"""
import openturns as ot
import numpy as np
import matplotlib
import pylab as plt
import warnings
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]


class PlotEE(object):

    """
    Plot elementary effects
    -----------------------

    The class perform the plot of sensitivity indices issued from Morris class.

    Parameters
    ----------
    morris: :class:`~otmorris.Morris`
        A Morris object.

    output_marginal: int
        Index of output marginal of interest.
        Default value is 0

    absolute_mean: bool
        Interest is mean of absolute elementary effects .
        Default value is True

    title: str
        Title for the graph

    """

    def __init__(self,
                 result,
                 output_marginal=0,
                 absolute_mean=True,
                 title="Elementary effects",
                 **kwargs):

        # set figure
        self._fig, self._ax = plt.subplots()
        # Check input object type
        if not (hasattr(result, 'getStandardDeviationElementaryEffects') and hasattr(result, 'getClassName')):
            raise TypeError(" `result` should be of class Morris ")
        if absolute_mean:
            mean = result.getMeanAbsoluteElementaryEffects(output_marginal)
        else:
            mean = result.getMeanElementaryEffects(output_marginal)
        sigma = result.getStandardDeviationElementaryEffects(output_marginal)
        dim = len(sigma)
        input_description = map(lambda x: "X" + str(x + 1), range(dim))
        # Plot effects
        self._ax.plot(mean, sigma, 'bo')
        # Annotate points
        dmu = (plt.np.max(mean) - plt.np.min(mean)) / len(mean)
        dsg = (plt.np.max(sigma) - plt.np.min(sigma)) / len(sigma)
        for i, txt in enumerate(input_description):
            self._ax.annotate(
                txt, (mean[i] + 0.05 * dmu, sigma[i] + 0.05 * dsg))
        self._ax.set_xlabel(r"$\boldsymbol{\mu}$", fontsize=14)
        self._ax.set_ylabel(r"$\boldsymbol{\sigma}$", fontsize=14)
        self._ax.grid(True)
        self._fig.suptitle(title, fontsize=18)

    def show(self, **kwargs):
        """
        Display the graph on screen.

        Parameters
        ----------
        kwargs:
            block: bool, optional
                If true (default), block until the graph is closed.

            These parameters are passed to matplotlib.pyplot.show()
        """
        self._fig.show(**kwargs)

    def save(self, fname, **kwargs):
        """
        Save the graph as file.

        Parameters
        ----------
        fname: bool, optional
            A string containing a path to a filename from which file format is deduced.

        kwargs:
            Refer to matplotlib.figure.Figure.savefig documentation for valid keyword arguments.
        """
        self._fig.savefig(fname, **kwargs)

    def getFigure(self):
        """
        Accessor to the underlying figure object.

        Refer to matplotlib.figure.Figure for further information.
        """
        return self._fig

    def getAxes(self):
        """
        Get the list of Axes objects.

        Refer to matplotlib.axes.Axes for further information.
        """
        return self._ax

    def close(self):
        """Close the figure."""
        plt.close(self._fig)
