import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from pathlib import Path

# Construct the paths
datasetName = "wine"
datasetPath = Path(__file__).parent.parent.parent / "datasets" / "csv" / datasetName
outputPath = Path(__file__).parent.parent.parent / "results" / datasetName

# Load the CSV file
data = pd.read_csv(str(datasetPath) + ".csv")

# Select only numeric columns
numeric_data = data.select_dtypes(include=["number"])

# Handle missing values (optional, depending on your data)
numeric_data = numeric_data.fillna(numeric_data.mean())

# Standardize the data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(numeric_data)

# Perform K-Means clustering
kmeans = KMeans(n_clusters=3, random_state=42)  # Set n_clusters to the desired number of clusters
kmeans.fit(scaled_data)

# Create a DataFrame for the cluster labels
cluster_column = pd.DataFrame(kmeans.labels_, columns=["Cluster"])

# Save only the cluster column to a new CSV file
cluster_column.to_csv(str(outputPath) + "_clusters.csv", index=False)

print("Clustering complete. Only cluster column saved to the CSV.")