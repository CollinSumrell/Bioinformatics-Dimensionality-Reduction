#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <filesystem>
#include <Eigen/Dense>

using namespace std;

typedef vector<vector<double>> Matrix;

Matrix centerData(const Matrix& data) {
    //this thang centers our data

    cout << "Centering data..." << "\n";
    int n = data.size(); //samples
    int m = data[0].size(); //features
    Matrix centered = data;
    
    for (int j = 0; j < m; j++) {
        //cout << j;

        double colSum = 0;
        for (int i = 0; i < n; i++) {
            colSum += data[i][j];
        }
        double mean = colSum / n;
        for (int i = 0; i < n; i++) {
            centered[i][j] -= mean;
        }
    }
    return centered;
}

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

//Use these later :) 

// Helper function to compute the magnitude of a vector
double vectorMagnitude(const vector<double>& vec) {
    double sumOfSquares = 0;
    for (double val : vec) {
        sumOfSquares += val * val;
    }
    return sqrt(sumOfSquares);
}

// Function to normalize a vector
vector<double> normalize(const vector<double>& vec) {
    double magnitude = vectorMagnitude(vec);
    vector<double> normalized(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        normalized[i] = vec[i] / magnitude;
    }
    return normalized;
}

// Function to perform PCA using the covariance matrix (simplified, assumes diagonalization is possible)
tuple<Matrix, vector<double>, Matrix> pca(const Matrix& data, int numComponents) {
    int n = data.size();  // Number of data points
    int m = data[0].size();  // Number of features

    Matrix centered = centerData(data);
    cout << "Data centered..." << "\n";

    Eigen::MatrixXd eigenData(n, m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            eigenData(i, j) = centered[i][j];
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
    Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors();

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

bool readCSV(const string& filename, Matrix& data) {
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

        bool firstCell = true; //flag to skip the first cell
        while (getline(ss, cell, ',')) { //split by comma
            if (firstCell) {
                firstCell = false; //skip the first cell in the row
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

int main() {
    // Example data matrix (4 samples, 3 features)
    // Matrix data = {
    //     {2.5, 2.4, 3.1},
    //     {0.5, 0.7, 2.2},
    //     {2.2, 2.9, 3.0},
    //     {1.9, 2.2, 2.8}
    // };


    Matrix data;
    if (!readCSV("datasets/csv/cells_CRC_0_A1_20220122143309lib1_z.csv", data)) {
        return 1;
    }

    cout << "data00: " << data[0][0] << "\n";
    cout << "Starting PCA" << "\n";

    

    int numComponents = 10;  // Number of principal components to keep
    // pca(data, numComponents);
    auto [eigenVectors, eigenValues, result] = pca(data, numComponents);

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

    return 0;
}