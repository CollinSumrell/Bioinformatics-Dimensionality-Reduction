### Useful Info: 

# barcodes.tsv is going to be our rownames, 
# features.tsv is going to be our col names

#matrix is our counts in the 1-indexed format of: 

# feature index (row) | barcode index (col) | count

### What this file does: 

# This 

import csv
from pathlib import Path

#construct the path
datasetPath = Path(__file__).parent.parent.parent / "datasets" 
binPath = Path(__file__).parent.parent.parent / "bin"

#loading barcodes (rows)
with open(datasetPath / "barcodes.tsv", "r") as f:
    print("Reading barcodes...")

    barcodes = [line.strip() for line in f]

    print("Finished reading barcodes.")

#loading features (cols)
with open(datasetPath / "features.tsv", "r") as f:
    print("Reading features...")

    features = []

    # Read each line and process it
    for line in f:
        parts = line.strip().split()  # Split the line into parts
        
        # Ensure there are at least 3 parts (to get the middle part)
        if len(parts) >= 3:
            middle_part = parts[1]  # The middle part is the second element
            features.append(middle_part)  # Add the middle part to the list

    print("Finished reading features")


nRow = len(barcodes)
nCol = len(features)

print(features)

print("nRow given by barcodes: " + str(nRow))
print("nCol given by features: " + str(nCol))

matrix = [[0] * nCol for _ in range(nRow)]

with open(datasetPath / "matrix.mtx", "r") as f:
    print("Reading matrix.mtx ...")

    while True:
        line = f.readline()
        if not line.startswith('%'):  # % lines are comments and are worthless, so we skip
            break
    
    #header matrix dimensions (features, samples, non-zero elements count)
    cols, rows, num_nonzeros = map(int, line.strip().split()) #strip removes any trailing spaces, split will split by tabs or whitespace
    
    if(nRow != rows or nCol != cols):
        print("COLUMN AND ROW LENGTH MISMATCH WHEN ASSEMBLING DATA, EXITING")
        exit()
    #read the file into the matrix
    for line in f:
        col, row, value = map(float, line.strip().split())
        matrix[int(row) - 1][int(col) - 1] = value  #convert to 0 indexing
    
    print("Finished reading matrix.mtx")

csvPath = binPath / "matrix.csv"

with open(csvPath, mode = "w", newline='') as csvfile:
    writer = csv.writer(csvfile)
    
    # Write the header (feature names as columns)
    writer.writerow(["Barcode"] + features)

    # Write each row (barcode) and its corresponding values
    for i, barcode in enumerate(barcodes):
        writer.writerow([barcode] + matrix[i])

    print(f"Matrix written to {binPath}")

