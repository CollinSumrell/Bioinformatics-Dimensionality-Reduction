import math 
import sys
import csv

from pathlib import Path

import time

# Add the utils folder to the Python path
#sys.path.append(str(Path(__file__).resolve().parent.parent / "utils"))

def transpose(matrix):
    #this exploits the way zip works by essentially turning each row
    #of the matrix into a backwards list, and then turning that into
    #its own list, which creates a transposed matrix. 

    return list(map(list,zip(*matrix)))

def mean(values): #simple mean function 
    return sum(values) / len(values)

#TODO: standardize data???? 

def centerData(data):
    #we center the data by subtracting the mean from everything 

    nRow = len(data) #samples
    nCol = len(data[0]) #features

    means = [mean([row[i] for row in data]) for i in range(nCol)]
    centeredData = [[data[i][j] - means[j] for j in range(nCol)] for i in range(nRow)]

    return centeredData,means #return our column meansm and our centeredData

def getCovarianceMatrix(data):
    #gets the covariance matrix , duh

    nRow = len(data) #samples
    nCol = len(data[0]) #features

    covMatrix = [[0 for _ in range(nCol)] for _ in range(nCol)] #this just initializes our array hence the "_"

    for i in range(nCol):
        for j in range(nCol):
            cov = sum(data[k][i] * data[k][j] for k in range(nRow)) / (nRow - 1)
            covMatrix[i][j] = cov

    #covMatrix = (transpose(data) @ data)/(nCol - 1)

    return covMatrix

#python already has a built in operator for this, but im going to use this code anyway so I 100% know whats happening

def matrixMultiply(A,B):
    result = [[0 for _ in range(len(B[0]))] for _ in range(len(A))]

    for i in range(len(A)):
        for j in range(len(B[0])):
            result[i][j] = sum(A[i][k] * B[k][j] for k in range(len(B)))
    return result

def vectorNorm(v):
    return math.sqrt(sum(x ** 2 for x in v))

def dotProduct(v1, v2):
    return sum(v1[i] * v2[i] for i in range(len(v1)))

def scalarMultiplyVector(scalar, vector):
    return [scalar * x for x in vector]

def vectorSubtract(v1, v2):
    return [v1[i] - v2[i] for i in range(len(v1))] #vectors must be same length

def powerIteration(matrix, numSimulations=100):

    n = len(matrix) #nRow, but in this case rows and columns are the same

    vectorBK = [1.0 for _ in range(n)] #vector of all ones
    for _ in range(numSimulations):
        # Calculate the matrix-by-vector product Ab
        vectorBK1 = [sum(matrix[i][j] * vectorBK[j] for j in range(n)) for i in range(n)]
        # Calculate the norm
        vectorBK1Norm = vectorNorm(vectorBK1)
        # Re-normalize the vector
        vectorBK = [x / vectorBK1Norm for x in vectorBK1]

        #this prevents the vector from growing too large

    # Compute the eigenvalue
    aB = [sum(matrix[i][j] * vectorBK[j] for j in range(n)) for i in range(n)]
    eigenValue = dotProduct(vectorBK, aB) / dotProduct(vectorBK, vectorBK)
    return eigenValue, vectorBK

def deflate(matrix, eigenValue, eigenVector):

    #removes an eigenValue and eigenVector pair from the matrix

    n = len(matrix)
    for i in range(n):
        for j in range(n):
            matrix[i][j] -= eigenValue * eigenVector[i] * eigenVector[j]
    return matrix

def readCSV(csvName, skipFirstCol):
    # Initialize an empty list to store rows
    data = []
    features = []
    barcodes = []

    # Open and read the CSV file
    with open(csvName, mode="r") as file:
        reader = csv.reader(file)

        features = next(reader)

        for row in reader:
            
            if skipFirstCol:
                barcodes.append(row[0])
                data.append(row[1:]) 
            else:
                data.append(row)  # Append each row as a list
    
    return data, barcodes, features

def pca(data, numComponents=None):
    #Center data before doing anything else
    centeredData, means = centerData(data)
    print("Data centered. Mean of each feature set to zero.")

    #Compute the covariance matrix
    covMatrix = getCovarianceMatrix(centeredData)
    
    print("Covariance matrix computed.")

    #Compute eigenvalues and eigenvectors using power iteration

    numFeatures = len(covMatrix)
    eigenValues = []
    eigenVectors = []

    covMatrixCopy = [row[:] for row in covMatrix]  # Make a copy for deflation

    for _ in range(numFeatures):
        val, vec = powerIteration(covMatrixCopy)
        eigenValues.append(val)
        eigenVectors.append(vec)
        covMatrixCopy = deflate(covMatrixCopy, val, vec)

    print("Eigenvalues and eigenvectors computed.")

    #Sort eigenvalues and eigenvectors in decreasing order

    eigenPairs = list(zip(eigenValues, eigenVectors))
    eigenPairs.sort(key=lambda x: x[0], reverse=True)
    eigenValues = [pair[0] for pair in eigenPairs]
    eigenVectors = [pair[1] for pair in eigenPairs]
    print("Eigenvalues and eigenvectors sorted in decreasing order.")

    #Select the top n_components if specified
    if numComponents is not None:
        eigenValues = eigenValues[:numComponents]
        eigenVectors = eigenVectors[:numComponents]

        print(f"Top {numComponents} principal components selected.")

    # Step 6: Project the data onto principal components
    projectedData = []
    for row in centeredData:
        projectedRow = [dotProduct(row, i) for i in eigenVectors]
        projectedData.append(projectedRow)
    print("Data projected onto principal components.")

    return projectedData, eigenValues, eigenVectors


if __name__ == "__main__":
    
    print("Starting PCA Python")

    # data = [
    #     [2.5, 2.4, 0.5],
    #     [0.5, 0.7, 0.3],
    #     [2.2, 2.9, 0.6],
    #     [1.9, 2.2, 0.4],
    #     [3.1, 3.0, 0.9]
    # ]

    printResults = False

    skipFirstCol = False
    datasetName = "wine.csv"
    datasetPath = Path(__file__).parent.parent.parent / "datasets" / "csv"
    binPath = Path(__file__).parent.parent.parent / "bin"


    data, barcodes, features = readCSV(datasetPath / datasetName, skipFirstCol)

    try: #converting our data to numeric since for some reason it comes in string form when we read it from the CSV
        data = [[float(value) for value in row] for row in data]
    except ValueError as e:
        print(f"Error converting data to numeric: {e}")
        sys.exit(1)

    start_time = time.perf_counter()

    # Perform PCA keeping all components
    projectedData, eigenValues, eigenVectors = pca(data,2)

    end_time = time.perf_counter()

    if(printResults):
        print("\nProjected Data:")
        for row in projectedData:
            print(row)

        print("\nEigenvalues:")
        print(eigenValues)

        print("\nEigenvectors:")
        for vec in eigenVectors:
            print(vec)

    elapsed_time = end_time - start_time
    print(f"Python took {elapsed_time:.6f} seconds to finish.")