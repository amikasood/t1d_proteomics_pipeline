import pandas as pd
import sqlite3
import numpy as np

def create_database(csv_path: str, db_path: str):
    '''
    Read the raw proteomics (csv) file provided in csv_path.
    Clean data
    Build database for further analysis

    Parameters
    ----------
    csv_path : path to csv input file
    db_path : Path to database file

    Returns
    -------
    None
    '''
    # Read csv proteomics file
    df = pd.read_csv(csv_path, header=1)
    #print(f"file path is {csv_path}")

    # Extract all annotations from the data
    annotations = df[['Protein ID', 'Gene', 'Annotated Gene', 'Annotated Matrisome Division', 'Annotated Matrisome Category']].copy()
    #print(annotations.head(5))

    # Extract expression data
    intensity_cols = [col for col in df.columns if 'MaxLFQ' in col]
    expression = df[['Protein ID'] + intensity_cols].copy()
    expression_long = expression.melt(id_vars=['Protein ID'], var_name='Sample', value_name='Intensity')

    # Clean data
    expression_long['Intensity'] = expression_long['Intensity'].replace(0, np.nan)
    expression_long['Sample ID'] = expression_long['Sample'].str.extract(r'(sample_\d+)')
    expression_long['Condition'] = np.where(expression_long['Sample ID'].isin([f"sample_{i:02d}" for i in [1, 2, 3, 5, 6, 7, 8, 9, 10]]), 'T1D', 'Control')

    # Create database
    conn = sqlite3.connect(db_path)
    annotations.to_sql('annotations', conn, if_exists='replace', index=False)
    expression_long[['Protein ID', 'Sample ID', 'Condition', 'Intensity']].to_sql('expression', conn, if_exists='replace', index=False) 

    print(f"Database created at {db_path}")