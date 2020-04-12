#tMFA

This module implements improved version of original thermodynamic meabolic flux analysis (tMFA).Christopher S. Henry, Linda J. Broadbelt, and Vassily Hatzimanikatis. "Thermodynamics-based metabolic flux analysis." Biophysical journal 92.5 (2007): 1792-1805. DOI: https://doi.org/10.1529/biophysj.106.093138

#Installation

Cloning the repository requires Git LFS to download some binary files. Git LFS can be found here (https://git-lfs.github.com/). To install, clone the repository using

```
git clone https://github.com/vishnumahamkali/tmfa.git
cd path/to/tmfa
git lfs install
git lfs pull
```
Please make sure to check if binary files in `Data` folder are properly downloaded for this to work. Please note, this module is developed and tested in Python 3.6.8 and should work on Python >= 3.6.

To install, simply go to the `tMFA` folder and

```
 python3 setup.py install
```
This module requires CoBRApy and it supports all COBRA compatiable solvers. It is recommended to use a commercial solver such as **GUROBI** or **CPLEX** to solve large MILP problems.

#Example script
To get started please see the script `example-tmfa.py`.

#Licence
The software in this repository is put under an APACHE-2.0 licensing scheme - please see the LICENSE file for more details.


