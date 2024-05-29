# The trait coding rule in phenotype space

05/28/24

Jianguo Wang (wangjg36@mail.sysu.edu.cn)
 

# UBHDD

UBHDD is a Python module for conducting the UBHDD algorithm, which is used for processing phenotype data frames and performing linear learning with uncorrelated features. It produces three output results: 
UBHDD_R2_R2SE_coef, containing R2, R2's standard error and coefficients of a trait in each row;
PG.csv
PNG.csv

## Installation

You can install the UBHDD module using pip:

```bash
pip install git+https://github.com/yourusername/UBHDD/UBHDD.git


## Usage
Once installed, you can import the UBHDD module in your Python scripts:

```python
from UBHDD import ubhdd

# Example usage
ubhdd.run_ubhdd('pheno_space.csv', output_folder="output", alpha=0.005, pcMatCutoff=None, n_processes=2)
