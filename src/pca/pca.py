##Some basics so I don't forget 




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

