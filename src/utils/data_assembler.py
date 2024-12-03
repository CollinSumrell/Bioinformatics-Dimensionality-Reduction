### Useful Info: 

# barcodes.tsv is going to be our rownames, 
# features.tsv is going to be our col names

#matrix is our counts in the 1-indexed format of: 

# feature index (row) | barcode index (col) | count

### What this file does: 

# This 

from pathlib import Path

#construct the path
datasetPath = Path(__file__).parent.parent.parent / "datasets" 

with open(datasetPath / "matrix.mtx", "r") as f:
    print("Reading matrix.mtx ...")

    while True:
        line = f.readline()
        if not line.startswith('%'):  # % lines are comments and are worthlesss
            break
    
    #header matrix dimensions (rows, columns, non-zero elements count)
    rows, cols, num_nonzeros = map(int, line.strip().split()) #strip removes any trailing spaces, split will split by tabs or whitespace
    
    matrix = [[0] * cols for _ in range(rows)] #big empty matrix
    
    #read the file into the matrix
    for line in f:
        row, col, value = map(float, line.strip().split())
        matrix[int(row) - 1][int(col) - 1] = value  # Convert to 0-based indexing


    print("Finished reading matrix.mtx")

