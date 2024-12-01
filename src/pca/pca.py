import math 

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


    covMatrix = [[0 for _ in range(nCol)] for _ in range(nRow)] #this just initializes our array hence the "_"

    for i in range(nRow):
        for j in range(nRow):
            cov = sum(data[k][i] * data[k][j] for k in range(nCol)) / (nCol - 1)
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
    return [v1[i] - v2[i] for i in range(len(v1))]

if __name__ == "__main__":
    
    #TODO: Load in a csv from 10x
    print("start")
    #hello