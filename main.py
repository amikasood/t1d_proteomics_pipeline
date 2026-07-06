import os
from src.data_handling import create_database

def main():
    # File paths
    csv_path = 'data/pancreas_proteomics.csv'
    db_path = 'data/pancreas_proteomics.db'

    if not os.path.exists(db_path):
        create_database(csv_path, db_path)
    else:
        print(f"Database already exists")

if __name__ == "__main__":
    main()