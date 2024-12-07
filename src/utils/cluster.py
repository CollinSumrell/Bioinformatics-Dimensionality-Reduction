import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

# Load the CSV file

from pathlib import Path

#construct the path

datasetName = "wine"

datasetPath = Path(__file__).parent.parent.parent / "datasets" / "csv" / datasetName

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

# Add the cluster labels to the original data
data["Cluster"] = kmeans.labels_

# Save the data with cluster labels to a new CSV file
data.to_csv(str(datasetPath) + "_clusters.csv", index=False)

print("Clustering complete. Results saved to 'clustered_data.csv'.")