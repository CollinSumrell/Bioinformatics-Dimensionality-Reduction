#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <filesystem>
#include <Eigen/Dense>
#include <chrono>

using namespace std;

typedef vector<vector<double>> Matrix;

// Function to calculate the covariance matrix
Matrix covarianceMatrix(const Matrix& data) {
    cout << "Finding covariance matrix...";

    int n = data.size();    // Number of samples
    int m = data[0].size(); // Number of features

    // Convert input data to Eigen matrix
    Eigen::MatrixXd eigenData(n, m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            eigenData(i, j) = data[i][j];
        }
    }

    // Compute the covariance matrix directly
    Eigen::MatrixXd covMatrix = (eigenData.transpose() * eigenData) / (n - 1);

    // Convert the result back to Matrix form out of Eigen form 
    Matrix result(m, std::vector<double>(m));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            result[i][j] = covMatrix(i, j);
        }
    }

    return result;
}

Matrix standardizeData(const Matrix& data) {
    int n = (int) data.size();
    int m = (int) data[0].size();

    Matrix standardized = data;

    for (int j = 0; j < m; ++j) {
        double colSum = 0.0;
        double colSumSq = 0.0;

        for (int i = 0; i < n; i++) {
            colSum += data[i][j];
            colSumSq += data[i][j] * data[i][j];
        }

        double mean = colSum / n;
        double var = (colSumSq / n) - (mean * mean);
        double stdDev = (var > 0) ? sqrt(var) : 1.0; // in case of zero variance

        for (int i = 0; i < n; i++) {
            standardized[i][j] = (data[i][j] - mean) / stdDev;
        }
    }

    return standardized;
}

// Function to perform PCA using the covariance matrix (simplified, assumes diagonalization is possible)
tuple<Matrix, vector<double>, Matrix> pca(const Matrix& data, int numComponents) {
    cout << "Start of PCA function" << endl;

    int n = data.size();  // Number of data points
    int m = data[0].size();  // Number of features

    Eigen::MatrixXd eigenData(n, m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            eigenData(i, j) = data[i][j];
        }
    }

    Eigen::MatrixXd covMatrix = (eigenData.transpose() * eigenData) / (n - 1);
    cout << "Covariance matrix calculated..." << "\n";

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(covMatrix);
    if (eigensolver.info() != Eigen::Success) {
        cerr << "Eigenvalue decomposition failed!" << endl;
        exit(1);
    }

    Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
    cout << "Computed eigenvalues..." << "\n";
    Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors();
    cout << "Computed eigenvectors..." << "\n";

    Eigen::VectorXd sortedEigenvalues = eigenvalues.reverse();
    Eigen::MatrixXd sortedEigenvectors = eigenvectors.rowwise().reverse();

    Eigen::MatrixXd topEigenvectors = sortedEigenvectors.leftCols(numComponents);

    Eigen::MatrixXd transformedData = eigenData * topEigenvectors;
    cout << "Projection completed..." << "\n";

    Matrix top_eigenvectors(m, vector<double>(numComponents));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < numComponents; ++j) {
            top_eigenvectors[i][j] = sortedEigenvectors(i, j);
        }
    }

    vector<double> eigenvalues_vec(numComponents);
    for (int i = 0; i < numComponents; ++i) {
        eigenvalues_vec[i] = sortedEigenvalues(i);
    }

    Matrix transformedDataMatrix(n, vector<double>(numComponents));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < numComponents; ++j) {
            transformedDataMatrix[i][j] = transformedData(i, j);
        }
    }

    return make_tuple(top_eigenvectors, eigenvalues_vec, transformedDataMatrix);
}
bool readCSV(const string& filename, Matrix& data, vector<string> barcodes, bool skipFirstCol) {
    cout << "Opening CSV: " << filename << "\n";

    ifstream file(filename); // Open the CSV file
    string line;

    if (!file.is_open()) {
        cerr << "Could not open the file!" << endl;
        return false;
    }

     if (getline(file, line)) {
        //skips first line, if we want to store feature names I can save here later
    }

    while (getline(file, line)) { 
        stringstream ss(line); //splits line via stringstream object
        string cell;
        vector<double> row; //stores each split in numeric


        bool skipFirstCol = true; //flag to skip the first cell
        while (getline(ss, cell, ',')) { //split by comma
            if (skipFirstCol) {

                barcodes.push_back(cell);

                skipFirstCol = false; //skip the first cell in the row
                continue;
            }
            try {
                double value = stod(cell);
                row.push_back(value);
            } catch (const invalid_argument& e) {
                cerr << "Invalid number: " << cell << endl;
                return false;
            }
        }

        if(!row.empty()){
            data.push_back(row);
        }
    }

    file.close(); // Close the file
    return true;
}

void writeToCSV(const Matrix& data, const string& filename) {
    ofstream file(filename);

    if (!file.is_open()) {
        cerr << "Failed to open the file: " << filename << endl;
        return;
    }

    // Add headers to the CSV file
    if (!data.empty()) {
        size_t numColumns = data[0].size();
        for (size_t i = 0; i < numColumns; ++i) {
            file << "PC" << (i + 1); // Generate headers like PC1, PC2, ...
            if (i < numColumns - 1) {
                file << ","; // Add a comma between headers, but not after the last one
            }
        }
        file << "\n"; // End the line after headers
    }

    // Write the data rows
    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i < row.size() - 1) {
                file << ","; // Add a comma between elements, but not after the last one
            }
        }
        file << "\n"; // End the line after each row
    }

    file.close();

    if (file.fail()) {
        cerr << "Error occurred while writing to the file: " << filename << endl;
    } else {
        cout << "Data successfully written to " << filename << endl;
    }
}

int main() {

    cout << "Starting PCA C++" << "\n";

    // Example data matrix (4 samples, 3 features)
    // Matrix data = {
    //     {2.5, 2.4, 3.1},
    //     {0.5, 0.7, 2.2},
    //     {2.2, 2.9, 3.0},
    //     {1.9, 2.2, 2.8}
    // };

    string datasetName = "wine";
    bool skipFirstCol = false; //if the dataset has built in row names, then skip 

    bool saveData = true; 
    bool printResults = false;

    int numComponents = 3;  // Number of principal components to keep

    vector<string> barcodes;
    Matrix data;

    if (!readCSV("datasets/csv/" + datasetName + ".csv", data, barcodes, skipFirstCol)) {
        return 1;
    }


    auto start = chrono::high_resolution_clock::now();

    data = standardizeData(data);

    // pca(data, numComponents);
    auto [eigenVectors, eigenValues, result] = pca(data, numComponents);

    auto end = chrono::high_resolution_clock::now();

    if(printResults == true){
        // Output the PCA result
        cout << "Eigenvectors:\n";
        for (const auto& vec : eigenVectors) {
            for (double val : vec) {
                cout << val << " ";
            }
            cout << endl;
        }

        cout << "\nEigenvalues:\n";
        for (double val : eigenValues) {
            cout << val << " ";
        }
        cout << endl;

        cout << "\nPCA Transformed Data (First " << numComponents << " Components):\n";
        for (const auto& row : result) {
            for (double value : row) {
                cout << value << " ";
            }
            cout << endl;
        }
    }

    string filename = "results/" + datasetName + "_" + std::to_string(numComponents) + "_PCs.csv";

    if(saveData){
        cout << "Saving data... \n";
        writeToCSV(result,filename);
    }

    std::chrono::duration<double> elapsed = end - start;

    // Output the elapsed time in seconds
    std::cout << "Function took " << elapsed.count() << " seconds to execute." << std::endl;

    return 0;
}