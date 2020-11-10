Thermodynamic Flux Variability Analysis & Sampling
==================================================

As demonstrated in the example script above, users can perform various types of
analyses including thermodynamic flux variability analysis (TVA) and flux
sampling. The workflow is described below,

.. math::

    \left( \mu \right)

Where :math:`\Sigma` is the cholesky matrix. The cholesky decomposition of a
positive-definite matrix is defined as

.. math::

    A = \Sigma \Sigma^{-1}

:math:`\Sigma` is a lower triangular matrix with real and positive diagonal
entries. If :math:`A` is positive semi-definite, then :math:`A` still has a
cholesky decomposition of the above form if the diagonal elements of cholesky
matrix are allowed to be zero. We use the algorithm described in [2]_ to compute
the nearest positive semi-definite covariance matrix.

We allow users to perform thermodynamic sampling in two ways: the box method and
sampling on the surface of an ellipsoid. In the box sampling method, we treat
formation energies as variables that are allowed to vary between mean and two
standard deviations as described in [3]_. This type of sampling over-estimates
the uncertainty of the reactions' Gibbs free energies. Thus, we
introduce another type of sampling method: sampling on the surface of an ellipsoid.
We achieve this by first sampling on the surface of the unit n-sphere and
transforming it to the surface of an n-ellipsoid. Sampling on the surface of the
unit n-sphere is described `here <https://mathworld.wolfram.com/HyperspherePointPicking.html>`_. 

We solve the tMFA problem for every sampled formation energy point on the
surface of the ellipsoid. Users can choose to sample until no futher improvement
is seen in reaction Gibbs free energies, after N sample points (user defined
cut-off), or we fit the user-defined number of samples to a `generalized extreme
value distribution
<https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.genextreme.html>`_
to predict the maximum and minimum reaction Gibbs free energies.

.. [2] Higham, Nicholas J.
    "Computing a Nearest Symmetric Positive Semidefinite Matrix."
    *Linear Algebra and Its Applications* 103 (May 1, 1988): 103–18.
    https://doi.org/10.1016/0024-3795(88)90223-6.

.. [3] Salvy, Pierre, Georgios Fengos, Meric Ataman, Thomas Pathier, Keng C Soh, and Vassily Hatzimanikatis.
    "PyTFA and MatTFA: A Python Package and a Matlab Toolbox for Thermodynamics-Based Flux Analysis."
    *Bioinformatics* 35, no. 1 (January 1, 2019): 167–69.
    https://doi.org/10.1093/bioinformatics/bty499.

