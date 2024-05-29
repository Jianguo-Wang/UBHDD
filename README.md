# The trait coding rule in phenotype space

05/28/24

Jianguo Wang (wangjg36@mail.sysu.edu.cn)
 

# UBHDD

UBHDD is a Python module for conducting the UBHDD (uncorrelation-based high-dimensional dependence) algorithm, which is designed for processing phenomic data and performing linear learning with uncorrelated features. It produces three output results:

1. `UBHDD_R2_R2SE_coef.csv`: Contains R-squared values, standard errors of R-squared, and coefficients for each trait in separate rows.
2. `PG.csv`: A file containing the genetic values for each trait.
3. `PNG.csv`: A file containing the non-genetic values for each trait.

## Installation

You can install these libraries using pip:

```bash
pip install git+https://github.com/Jianguo-Wang/UBHDD.git

```

## Usage

Once installed, you can import the UBHDD module in your Python scripts and use the `run_ubhdd` function:

```python
from UBHDD import ubhdd

# Example usage
ubhdd.run_ubhdd('pheno_space.csv', output_folder="output", alpha=0.005, pcMatCutoff=None, n_processes=2)
```

The `run_ubhdd` function takes the following arguments:

- `'pheno_space.csv'` is the path to your phenotype data file in CSV format.
- `output_folder="output"` specifies the folder where the output files will be saved.
- `alpha=0.005` sets the significance level for determining the uncorrelation threshold.
- `pcMatCutoff=None` set the uncorrelation threshold directly.
- `n_processes=2` specifies the number of parallel processes to use for computation.

After running the `run_ubhdd` function, the following output files will be generated in the specified `output` folder:

1. `UBHDD_R2_R2SE_coef.csv`: Contains the R-squared values, standard errors of R-squared, and coefficients for each trait.
2. `PG.csv`: A file containing the genetic values for each trait.
3. `PNG.csv`: A file containing the non-genetic values for each trait.


```
