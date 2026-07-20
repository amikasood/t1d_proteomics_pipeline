import pandas as pd
import numpy as np
import sqlite3
from scipy import stats
import statsmodels.stats.multitest as mt

def run_differential_expression(db_path):
    '''
    Run differential expression analysis (MNAR data)

    Parameters
    ----------
    db_path : path to the database file

    Returns
    -------
    DEP : A DataFrame of differentially expressed proteins
    results_df : A Dataframe of Protein ID, Log2FC, p-value and adjusted p-value 
    '''
    with sqlite3.connect(db_path) as conn:
        df = pd.read_sql("SELECT * FROM imputed_MNAR", conn)
    
    df_pivot = df.pivot(index='Protein ID', columns='Sample ID', values='Intensity')

    results = []
    for protein in df_pivot.index:
        index = df_pivot.index.get_loc(protein)
        t1d = df_pivot.iloc[index, list(range(9))].dropna()
        ctrl = df_pivot.iloc[index, list(range(9, 18))].dropna()

        _, p_val = stats.ttest_ind(t1d, ctrl, equal_var=False)
        fc = t1d.mean() - ctrl.mean()

        results.append((protein, fc, p_val))

    results_df = pd.DataFrame(results, columns=['Protein ID', 'Log2FC', 'p_value'])
    _, results_df['adj_p_value'], _, _ = mt.multipletests(results_df['p_value'], method='fdr_bh')

    dep = results_df[(results_df['adj_p_value'] < 0.05) & (results_df['Log2FC'].abs() > 0.58)]

    return dep, results_df
