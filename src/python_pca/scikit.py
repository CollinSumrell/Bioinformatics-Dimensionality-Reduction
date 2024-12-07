import pandas as pd
import numpy as np
import time
from sklearn.decomposition import PCA as SklearnPCA
from sklearn.preprocessing import StandardScaler
from matplotlib import pyplot as plt
from prince import PCA as PrincePCA
from pathlib import Path

datasetName = "wine.csv"
datasetPath = Path(__file__).parent.parent.parent / "datasets" / "csv"

numComponents = 3  # This should be used with PCA, not StandardScaler

# Load the dataset from CSV
file_path = datasetPath / datasetName  # Replace with your dataset path
data = pd.read_csv(file_path)

# Ensure numeric data only
data = data.select_dtypes(include=[np.number])

# Standardize the dataset (do not pass numComponents to StandardScaler)
scaler = StandardScaler()
data_scaled = scaler.fit_transform(data)

# Function to measure PCA execution time and return transformed data
def measure_pca_time(pca_library, data, library_type="sklearn"):
    start_time = time.time()
    if library_type == "sklearn":
        reduced_data = pca_library.fit_transform(data)
    elif library_type == "prince":
        reduced_data = pca_library.fit(data).row_coordinates(data)
    elapsed_time = time.time() - start_time
    return reduced_data, elapsed_time

# PCA using sklearn
sklearn_pca = SklearnPCA(n_components=numComponents)  # Correctly pass numComponents to PCA
sklearn_pca_data, sklearn_time = measure_pca_time(sklearn_pca, data_scaled)

# PCA using prince
prince_pca = PrincePCA(n_components=numComponents)  # Correctly pass numComponents to PCA
prince_pca_data, prince_time = measure_pca_time(prince_pca, pd.DataFrame(data_scaled), library_type="prince")

# Print execution times
print(f"Sklearn PCA Time: {sklearn_time:.4f} seconds")
print(f"Prince PCA Time: {prince_time:.4f} seconds")