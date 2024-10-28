User manual
===========
In the context of sensitivity analysis, we often have to deal with high dimensional inputs and heavy CPU time codes.
Thus usual sensitivity techniques based on ANOVA reach their limits as they require to discretize conditional variance, which needs several samples
to get an accurate approximation.

To tackle this kind of problems, screening methods could be applied in order to get a qualitative sensitivity estimate.
The Morris method address the needs.

The method focuses on the notion of elementary effects and is known to require very few simulations to get an `accurate estimate`
of the influent factors.

Roughly speeking, the method relies on One At Time designs (OAT) and acts as follows:

- The input design space is discretized (in a p-levels grid), of step :math:`\Delta`;
- We randomly choose a starting `point` in this grid;
- We randomly select a direction and thus we get the new point,
- We iterate the previous process on the :math:`p-1` remaining directions to get a full `path`, where :math:`p` is the input dimension.
  Note that :math:`p+1` experiments are needed to define this `path`. (See hereafter an example in case :math:`p=2`)

.. image:: scheme-path-morris.png
   :height: 350 px
   :width: 600 px
   :align: center

- From this path (:math:`X_p`), we compute the response answer :math:`Y_p=f(X_p)`;
- It is easy to see that the difference between two consecutive points of this `path` represents the elementary effect
  relative to the choosen direction.
  Indeed we compute both :math:`dX_p` and :math:`dY_p` where :math:`dX_p` represents the difference between
  two consecutive elements of :math:`X_p`.
  We deduce elementary effects from these vectors of size :math:`p` by solving the linear system
  :math:`X_p\ ee = Y_p` (:math:`ee` are the elementary effects)
- We iterate the steps 2-5 :math:`r` times in order to get `r` replicates of the elementary effects.
  Here after an illustration  in case :math:`p=2, r=5`.

.. image:: scheme-path-examples.png
   :height: 350 px
   :width: 400 px
   :align: center

If we note :math:`ee_i^k` the k-th computed elementary effects associated to the i-th input marginal, it follows that from the r-sample of elementary effects,
we get :math:`\mu_i^*, \sigma_i` respectively the absolute mean and the standard deviation of the elementary effects:


.. math::
   \mu_i^* = \frac{1}{r} \sum_{k=1}^{r} |ee_i^k|

   \sigma_i = \sqrt{ \frac{1}{r} \left(\sum_{k=1}^{r} ee_i^k - \mu_i \right)^2 }

These are the measure used to get a qualitative approch of the sensitivity. In the original Morris implementation, the mean of elementary effects :math:`\mu_i = \frac{1}{r} \sum_{k=1}^{r} ee_i^k` was used but it lacks of precision due to some potential sign changes. The previous values could be interpreted as follows:
  - :math:`\mu_i^*` measures sensitivity. Important values highlight important effects and thus that model is sensitive to input variations,
  - :math:`\sigma_i` measures the interactions and non linearity effects. Important values could explains non linear effects or interactions but it is impossible to make the distinction between the two cases.


In engineering application litterature, other interpretations could be found based on the quantity :math:`\rho_i = \frac{\sigma_i}{\mu_i^*}`:
  - If :math:`\rho_i \leq 0.1` the i-th variable has almost linear effects,
  - If :math:`0.1 \leq \rho_i \leq 0.5` the i-th variable has monotonic effects,
  - If :math:`0.5 \leq \rho_i \leq 1` the i-th variable has quasi-monotonic effects,
  - If :math:`\rho_i \geq 1` the i-th variable has non-linear and non-monotonic effects


To conclude, this module allows to estimate the previous sensitivity measures (both :math:`\mu, \mu^*, \sigma`) starting both from a `p-level` grid or an `LHS` experiment. It allows also to get response model outside the library and finally plot the sensitivity to get a qualitative estimate.


Reference
---------
 - Campolongo, F., S. Tarantola and A. Saltelli. (1999). "Tackling quantitatively
   large dimensionality problems.". Computer Physics Communication 1999: 75–85.
   `pdf <http://www.andreasaltelli.eu/file/repository/Campolongo_Saltelli1_s2_0_S0010465598001659_main.pdf>`__
 - Morris, M.D. (1991). "Factorial Sampling Plans for Preliminary Computational
   Experiments" (PDF). Technometrics 33: 161–174.
   `pdf <http://abe.ufl.edu/Faculty/jjones/ABE_5646/2010/Morris.1991%20SA%20paper.pdf>`__
 - Campolongo, Cariboni, Saltelli, F., J.and A. (2003). "Sensitivity analysis: the
   Morris method versus the variance based measures. Submitted to Technometrics"
   (PDF).
 - Andrea Saltelli, Stefano Tarantola,Francesca Campolongo and Marco Ratto
   (2004). "Sensitivity analysis in practice a guide to assessing scientific
   models". John Willy & sons
   `pdf <http://www.andreasaltelli.eu/file/repository/SALTELLI_2004_Sensitivity_Analysis_in_Practice.pdf>`__


Experiments for Morris
----------------------

.. currentmodule:: otmorris
.. autosummary::
    :toctree: _generated/
    :template: class.rst_t

    MorrisExperiment
    MorrisExperimentGrid
    MorrisExperimentLHS


Morris screening method
-----------------------
.. currentmodule:: otmorris
.. autosummary::
    :toctree: _generated/
    :template: class.rst_t

    Morris


Morris function
---------------
.. currentmodule:: otmorris
.. autosummary::
    :toctree: _generated/
    :template: class.rst_t

    MorrisFunction
