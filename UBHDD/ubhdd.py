import multiprocessing
import os
import numpy as np
import pandas as pd
from glmnet import ElasticNet
import scipy.stats as stats

def compute_correlation_matrix(pheno_space):
    """Compute the correlation matrix."""
    pheno_space3 = pheno_space.iloc[:, 1:].to_numpy()
    return np.corrcoef(pheno_space3, rowvar=False)

def linear_learn_unCor(args):
    """Function for linear learning using uncorrelated features."""
    i, pheno_space3, pheno_space_corMat3, pcMatCutoff = args
    m = ElasticNet(n_splits=10, n_lambda=200, scoring='r2')
    focalUnCorInd = np.nonzero(np.abs(pheno_space_corMat3[:, i]) < pcMatCutoff)[0]
    fit = m.fit(pheno_space3[:, focalUnCorInd], pheno_space3[:, i])
    lambda_max_ind = np.argwhere(fit.lambda_path_ == fit.lambda_max_)[0, 0]
    coef = np.zeros(pheno_space3.shape[1])
    coef[focalUnCorInd] = fit.coef_path_[:, lambda_max_ind]
    return [fit.cv_mean_score_[lambda_max_ind], fit.cv_standard_error_[lambda_max_ind],
            fit.intercept_path_[lambda_max_ind]] + coef.reshape(1, pheno_space3.shape[1]).tolist()[0]

def run_ubhdd(file_path, output_folder="output", alpha=0.005, pcMatCutoff=None, n_processes=80):
    """Run UBHDD algorithm."""
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # Load phenotype data
    pheno_space = pd.read_csv(file_path)
    pheno_space_ColName = pheno_space.columns[1:]
    pheno_space_RowName = pheno_space.iloc[:, 0]
    
    # Compute the correlation matrix
    pheno_space_corMat3 = compute_correlation_matrix(pheno_space)

    # Define parameters
    n_sample, n_trait = pheno_space.shape[0], pheno_space.shape[1] - 1

    # Calculate the critical value using the t-distribution
    if pcMatCutoff is None:
        critical_value = stats.t.ppf(1 - alpha / (n_trait - 1), n_sample - 2)
        pcMatCutoff = np.sqrt(critical_value ** 2 / ((n_sample - 2) + critical_value ** 2))

    # Use multiprocessing to parallelize the computation
    args_list = [(i, pheno_space3, pheno_space_corMat3, pcMatCutoff) for i in range(n_trait)]
    with multiprocessing.Pool(processes=n_processes) as pool:
        result = pool.map(linear_learn_unCor, args_list)

    # Save the results to a DataFrame
    UBHDD_R2_R2SE_coef = pd.DataFrame(result, columns=['R2', 'R2SE', 'Intercept'] + list(pheno_space_ColName),
                                       index=pheno_space_ColName)
    UBHDD_R2_R2SE_coef.to_csv(os.path.join(output_folder, 'UBHDD_R2_R2SE_coef.csv'))

    # Compute predicted genetic values and residuals
    PG = pheno_space3 @ UBHDD_R2_R2SE_coef.iloc[:, 2:].to_numpy().T
    PNG = pheno_space3 - PG

    # Save the results to CSV files
    PG_df = pd.DataFrame(PG, columns=list(pheno_space_ColName), index=pheno_space_RowName)
    PNG_df = pd.DataFrame(PNG, columns=list(pheno_space_ColName), index=pheno_space_RowName)

    PG_df.to_csv(os.path.join(output_folder, 'PG.csv'))
    PNG_df.to_csv(os.path.join(output_folder, 'PNG.csv'))

