from pathlib import Path
import pandas as pd

def append_column_to_csv(source_csv, target_csv):
    # Read the source and target CSVs
    source_data = pd.read_csv(source_csv)
    target_data = pd.read_csv(target_csv)
    
    # Ensure the source CSV has exactly one column
    if source_data.shape[1] != 1:
        raise ValueError("Source CSV must have exactly one column.")
    
    # Check if the last column of the target CSV is the same as the source column
    if target_data.columns[-1] == "Cluster":
        print("The target CSV already contains a 'Cluster' column as the last column. No changes made.")
        return
    
    # Append the column from source to target
    source_data.columns = ["Cluster"]  # Ensure the column is named "Cluster"
    target_data = pd.concat([target_data, source_data], axis=1)
    
    # Overwrite the target CSV
    target_data.to_csv(target_csv, index=False)
    print(f"Appended column from {source_csv} to {target_csv} and saved.")

datasetName = "proteindata"

sourcePath = Path(__file__).parent.parent.parent / "results" / f"{datasetName}_clusters.csv"
targetPath = Path(__file__).parent.parent.parent / "results" / f"{datasetName}_3_PCs.csv"
append_column_to_csv(sourcePath, targetPath)