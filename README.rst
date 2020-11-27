========
multiTFA
========

.. image:: https://img.shields.io/pypi/v/multitfa.svg
   :target: https://pypi.org/project/multitfa/
   :alt: Current PyPI Version

.. image:: https://img.shields.io/pypi/pyversions/multitfa.svg
   :target: https://pypi.org/project/multitfa/
   :alt: Supported Python Versions

.. image:: https://img.shields.io/pypi/l/multitfa.svg
   :target: https://www.apache.org/licenses/LICENSE-2.0
   :alt: Apache Software License Version 2.0

.. image:: https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg
   :target: .github/CODE_OF_CONDUCT.md
   :alt: Code of Conduct

.. image:: https://github.com/biosustain/multitfa/workflows/CI-CD/badge.svg
   :target: https://github.com/biosustain/multitfa/workflows/CI-CD
   :alt: GitHub Actions

.. image:: https://codecov.io/gh/biosustain/multitfa/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/biosustain/multitfa
   :alt: Codecov

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/ambv/black
   :alt: Code Style Black

.. image:: https://readthedocs.org/projects/multitfa/badge/?version=latest
   :target: https://multitfa.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. summary-start

We present multiTFA, a multivariate thermodynamics-based metabolic flux analysis
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

    pip install multitfa

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
