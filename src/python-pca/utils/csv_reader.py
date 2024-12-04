import csv
from pathlib import Path

#construct the path
datasetPath = Path(__file__).parent.parent.parent / "datasets" 
binPath = Path(__file__).parent.parent.parent / "bin"

def readCSV(csvName):
    # Initialize an empty list to store rows
    data = []
    features = []
    barcodes = []

    # Open and read the CSV file
    with open(binPath / csvName, mode="r") as file:
        reader = csv.reader(file)

        features = next(reader)

        for row in reader:
            barcodes.append(row[0])
            data.append(row[1:])  # Append each row as a list
    
    return data, barcodes, features