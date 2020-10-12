=================
Multivariate tMFA
=================

.. image:: https://img.shields.io/pypi/v/tmfa.svg
   :target: https://pypi.org/project/tmfa/
   :alt: Current PyPI Version

.. image:: https://img.shields.io/pypi/pyversions/tmfa.svg
   :target: https://pypi.org/project/tmfa/
   :alt: Supported Python Versions

.. image:: https://img.shields.io/pypi/l/tmfa.svg
   :target: https://www.apache.org/licenses/LICENSE-2.0
   :alt: Apache Software License Version 2.0

.. image:: https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg
   :target: .github/CODE_OF_CONDUCT.md
   :alt: Code of Conduct

.. image:: https://github.com/biosustain/multivariate-tmfa/workflows/CI-CD/badge.svg
   :target: https://github.com/biosustain/multivariate-tmfa/workflows/CI-CD
   :alt: GitHub Actions

.. image:: https://codecov.io/gh/biosustain/multivariate-tmfa/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/biosustain/multivariate-tmfa
   :alt: Codecov

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/ambv/black
   :alt: Code Style Black

.. image:: https://readthedocs.org/projects/tmfa/badge/?version=latest
   :target: https://tmfa.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. summary-start

We present ``tmfa``, a multivariate thermodynamics-based metabolic flux analysis
package for Python. The framework takes advantage of the reactions' Gibbs free
energy covariance matrix to tightly constrain metabolic models using
thermodynamic constraints. It represents an improvement over a previous
thermodynamic metabolic flux analysis (tMFA) method described in [1]_.

This implementation requires a COBRA model as input, as well as additional
compartment conditions and metabolite concentrations. It allows user to perform
various thermodynamic analyses on COBRA models, such as thermodynamic metabolic
flux analysis, variability analysis, or flux sampling.  Please see below for
further details.

.. [1] Henry, Christopher S., Linda J. Broadbelt, and Vassily Hatzimanikatis.
    "Thermodynamics-Based Metabolic Flux Analysis."
    *Biophysical Journal* 92, no. 5 (March 1, 2007): 1792–1805.
    https://doi.org/10.1529/biophysj.106.093138.

Install
=======

It's as simple as:

.. code-block:: console

    pip install tmfa

We highly recommend the installation of a commercial mathematical optimization
solver, like `GUROBI <https://www.gurobi.com/>`_ or `CPLEX
<https://www.ibm.com/analytics/cplex-optimizer>`_.

Usage
=====

Please take a look at the Python scripts and Jupyter notebooks in the
``examples`` directory.

Copyright
=========

* Copyright © 2018, Novo Nordisk Foundation Center for Biosustainability.
* Free software distributed under the `Apache Software License 2.0
  <https://www.apache.org/licenses/LICENSE-2.0>`_.

.. summary-end

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

