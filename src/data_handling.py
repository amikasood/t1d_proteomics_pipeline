import pandas as pd
import sqlite3
import numpy as np
import mygene

def create_database(csv_path: str, db_path: str):
    '''
    Read the raw proteomics (csv) file provided in csv_path.
    Clean data
    Add annotations
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

    # Add annotations 
    # Add panther annothations to the non-matrisome proteins
    non_matrisome_annotation = pd.read_csv("../data/pantherGeneList.txt", sep="\t", header=False)
    df['Annotated Matrisome Category'] = df['Annotated Gene'].map(non_matrisome_annotation.set_index('Annotated Gene')['Annotated Matrisome Category']).fillna(df['Annotated Matrisome Category']).fillna("Uncalssified")

    # Add GO:CC (Cellular component annotations)
    mg = mygene.MyGeneInfo()
    go_cc_ann = mg.querymany(df['Annotated Gene'].to_list(), scopes='uniprot', fields='go', species='human', as_dataframe=True)
    go_cc_ann = go_cc_ann.reset_index()
    go_cc_ann_clean = go_cc_ann.drop_duplicates(subset='query', keep='first')
    go_cc_ann_clean = go_cc_ann_clean.reset_index(drop=True)
    go_cc_ann_clean['GO CC'] = go_cc_ann_clean['go.CC'].apply(categorize_go_cc)
    
    # Extract all annotations from the data
    annotations = df[['Protein ID', 'Gene', 'Annotated Gene', 'Annotated Matrisome Division', 'Annotated Matrisome Category']].copy()
    annotations['GO CC'] = go_cc_ann_clean['GO CC'].values
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

# Categorization function using GO CC
def categorize_go_cc(cc_data):
    # Check for NaN or None
    if isinstance(cc_data, float) or cc_data is None:
        return 'Unclassified'
    # Check for empty list
    if isinstance(cc_data, list) and len(cc_data)==0:
        return 'Unclassified'
    # Single dictionaries
    if isinstance(cc_data, dict):
        cc_data = [cc_data]
    
    term_names = " ".join([item.get('term', '').lower() for item in cc_data if isinstance(item, dict)])

    if 'mitochondri' in term_names:
        return 'Mitochondrial Proteins'
    elif 'nucleus' in term_names or 'nuclear' in term_names:
        return 'Nuclear Proteins'
    elif 'endoplasmic reticulum' in term_names or 'golgi' in term_names:
        return 'ER/Golgi Proteins'
    elif 'cytoskelet' in term_names or 'actin' in term_names or 'microtubule' in term_names:
        return 'Cytoskeletal Proteins'
    elif 'membrane' in term_names:
        return 'Membrane Proteins'
    elif 'cytoplasm' in term_names or 'cytosol' in term_names:
        return 'Cytosolic Proteins'
    elif 'secretory granule' in term_names or 'vesicle' in term_names or 'exosome' in term_names:
        return 'Secretory/Vesicular Proteins'
    elif 'extracellular' in term_names or 'secretes' in term_names:
        return 'ECM/secreted proteins'
    elif 'lysosom' in term_names or 'endosom' in term_names or 'autophagosom' in term_names:
        return 'Lysosomal/Endosomal Proteins'
    elif 'ribosom' in term_names or 'proteasom' in term_names:
        return 'Ribosomal/Proteasomal Proteins'
    elif 'plasma membrane' in term_names or 'cell surface' in term_names:
        return 'Plasma membrane Proteins'
    else:
        return 'Other'