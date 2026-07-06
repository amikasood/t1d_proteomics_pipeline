import pandas as pd
import numpy as np
import sqlite3
from sklearn.impute import KNNImputer

def load_and_pivot(db_path):
    '''
    Read data from database and reshape it from a long format to a wide format 
    Create a lookup key matching Sample ID and Condition

    Parameters
    ----------
    db_path :  Path to the database

    Returns
    -------
    df_pivot : DataFrame of Intensity values
    conditions_map : Dictionary that maps Sample ID to Condition
    '''
    with sqlite3.connect(db_path) as conn:
        df = pd.read_sql("SELECT [Protein ID], [Sample ID], Condition, Intensity FROM expression", conn)
    
    df_pivot = df.pivot(index='Protein ID', columns='Sample ID', values='Intensity')

    conditions_df = df[['Sample ID', 'Condition']].drop_duplicates()
    condition_map = dict(zip(conditions_df['Sample ID', 'Condition']))

    return df_pivot, condition_map

def filter_missing_values(df, condition_map, threshold):
    '''
    Filter proteins with missingness more than threshold

    Parameters
    ----------
    df : DataFrame of intensity values
    condition_map : Dictionary that maps Sample ID to Condition
    threshold : maximum allowed missingness

    Returns
    -------
    df_filter : DataFrame of intensity values. Proteins filtered based on threshold
    '''
    t1d_id = [id for id, cond in condition_map.items() if cond == 'T1D']
    ctrl_id = [id for id, cond in condition_map.items() if cond == 'Control']

    df['T1D_count'] = df[t1d_id].notna().sum(axis=1)
    df['Ctrl_count'] = df[ctrl_id].notna().sum(axis=1)

    df_filter = df[(df['T1D_count'] >= threshold) | (df['Ctrl_count'] >= threshold)]
    df_filter = df_filter.drop(columns=['T1D_count', 'Ctrl_count'])

    return df_filter

def normalize(df):
    '''
    Normalize and center data

    Parameters
    ----------
    df : DataFrame of intensity values

    Returns
    -------
    df_norm : normalized and centered DataFrame of intensity values
    '''
    df_log = np.log2(df)
    sample_median = np.median(df_log, axis=0)
    global_median = sample_median.median()
    df_norm = df_log.subtract(sample_median, axis=1) + global_median
    return df_norm

def imputation(df, n_neighbors):
    '''
    Impute missing data with kNN inputation

    Parameters
    ----------
    df : DataFrame of intensity values
    n_neighbors : number of similar rows the algorithm will examine

    Returns
    -------
    df_imputed : Imputed DataFrame of intensity values
    '''
    impute = KNNImputer(n_neighbors=n_neighbors, weights='distance')
    imputed_array = impute.fit_transform(df.T)
    df_imputed = pd.DataFrame(imputed_array, index=df.T.index, columns=df.T.columns).T
    return df_imputed

def pre_processing(db_path):
    print("Starting data processing and imputation.... ")
    df_pivot, condition_map = load_and_pivot(db_path)
    df_filter = filter_missing_values(df_pivot, condition_map, 5)
    df_norm = normalize(df_filter)

    df_imputed_KNN = imputation(df_norm, n_neighbors=3)

    shifted_min = df_norm.min().min() - 1.0
    df_imputed_mnar =  df_norm.fillna(shifted_min)

    with sqlite3.connect(db_path) as conn:
        long_KNN = df_imputed_KNN.reset_index().melt(id_vars='Protein ID', var_name='Sample ID', value_name='Intensity')
        long_KNN.to_sql('imputed_KNN', conn, if_exists='replace', index=False)
        long_mnar = df_imputed_mnar.reset_index().melt(id_vars='Protein ID', var_name='Sample ID', value_name='Intensity')
        long_mnar.to_sql('imputed_MNAR', conn, if_exists='replace', index=False)
    
    print("Data saved to ../pancreas_proteomics.db")
