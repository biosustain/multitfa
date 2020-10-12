# multitFA

This module implements improved version of original thermodynamic meabolic flux
analysis (tMFA).Christopher S. Henry, Linda J. Broadbelt, and Vassily
Hatzimanikatis. "Thermodynamics-based metabolic flux analysis." Biophysical
journal 92.5 (2007): 1792-1805. DOI:
https://doi.org/10.1529/biophysj.106.093138.

We present multiTFA, framework that takes advantage of the covariance matrix to
tightly constrain the metabolic models using thermodynamic constraints. This
implementation requires COBRA model and compartment and metabolite information
as input. This framework allows user to perform various thermodynamic analyses
on COBRA models including thermodynamic metabolic flux analysis, variability
analysis, sampling (Please see below for further details). 

## Installation

Cloning the repository requires Git LFS to download some binary files. Git LFS
can be found [here](https://git-lfs.github.com/). To install, clone the
repository using

```
git clone https://github.com/vishnumahamkali/tmfa.git
cd path/to/tmfa
git lfs install
git lfs pull
```
Please make sure to check if binary files in `Data` folder are properly
downloaded for this to work. Please note, this module is developed and tested in
Python 3.6.8 and should work on Python >= 3.6.

To install, simply go to the `tMFA` folder and

```
 python3 setup.py install
```
This module requires COBRApy and it supports all COBRA compatible solvers. It is
recommended to use a commercial solver such as **GUROBI** or **CPLEX** to solve
large MILP problems.

## Example script

To get started please see the script `example-tmfa.py`.

## Licence

The software in this repository is put under an APACHE-2.0 licensing scheme - please see the LICENSE file for more details.

## Thermodynamic variability analysis (TVA) & sampling

As demonstrated in the example script above, users can perform various types of analyses including TVA and sampling. The workflow is described below,


$(𝜇^−𝜇)<sup>T</sup> 𝛴<sup>−1</sup> (𝜇^−𝜇) ≤ \chi<sup>2</sup><sub>n,95%</sub>$

Where 𝛴 is cholesky matrix. The cholesky decomposition of a positive-definite matrix is defined as


$A = 𝛴 𝛴<sup>-1</sup>$

𝛴 is a lower triangular matrix with real and positive diagonal entries. If A is positive semi-definite, then A still has cholesky decomposition of above form if the diagonal elements of cholesky matrix is allowed to be zero. We use algorithm described by [Higham et.al](https://doi.org/10.1016/0024-3795(88)90223-6) to compute the nearest positive semi-definite covariance matrix. 

We allow users to perform thermodynamic sampling in two ways, the box method and sampling on the surface of ellipsoid. In the box sampling method, we treat formation energies as variables that are allowed to vary between mean and two standard deviations as described by [Salvy et.al](https://doi.org/10.1093/bioinformatics/bty499). This type of sampling is not ideal (Please refer to our manuscript for futher details). Thus, we introduce another type of sampling method, sampling on the surface of ellipsoid. We achieve this by first sampling on the surface of the unit n-sphere and transforming it to the surface of the n-ellipsoid. Sampling on the surface of n-sphere is described [here](https://mathworld.wolfram.com/HyperspherePointPicking.html). 

We solve tMFA problem for every sampled formation energy point on surface of ellipsoid. User can choose to sample until no futher improvement is seen in reaction Gibbs free energies after N sample points (user defined cut-off) or we fit the user-defined number of samples to [generalized extreme value distribution](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.genextreme.html) to predict the maximum and minimum reaction Gibbs free energies. 

## Debugging

When the model status is 'infeasible', users can find the constraints that
render the model infeasible. Please note this functionality works only when
using [GUROBI](https://www.gurobi.com) solver. An example of how to use this
functionality is shown below.

```
while np.isnan(model.slim_optimize()):
	model.solver.problem.computeIIS()
	for c in model.solver.problem.getConstrs():
		if c.IISConstr:
            print(c)
```

Please note there can be various combination of constraints that cause problems,
users are advised to account for biological meaning when removing a constraint. 
