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
    int n = data.size();
    int m = data[0].size();
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
Matrix covariance_matrix(const Matrix& data) {
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

    // Convert the result back to std::vector<std::vector<double>>
    Matrix result(m, std::vector<double>(m));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            result[i][j] = covMatrix(i, j);
        }
    }

    return result;
}

// Helper function to compute the magnitude of a vector
double vector_magnitude(const vector<double>& vec) {
    double sum_of_squares = 0;
    for (double val : vec) {
        sum_of_squares += val * val;
    }
    return sqrt(sum_of_squares);
}

// Function to normalize a vector
vector<double> normalize(const vector<double>& vec) {
    double magnitude = vector_magnitude(vec);
    vector<double> normalized(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        normalized[i] = vec[i] / magnitude;
    }
    return normalized;
}

// Function to perform PCA using the covariance matrix (simplified, assumes diagonalization is possible)
tuple<Matrix, vector<double>, Matrix> pca(const Matrix& data, int num_components) {
    int n = data.size();  // Number of data points
    int m = data[0].size();  // Number of features

    cout << "actually starting pca for real this time (maybe)" << "\n";

    // Step 1: Center the data
    Matrix centered = centerData(data);

    cout << "Data centered..." << "\n";

    // Step 2: Compute covariance matrix
    Matrix cov_matrix = covariance_matrix(centered);

    cout << "Covariance matrix computed..." << "\n";

    // Step 3: Eigenvalue decomposition (simplified, uses basic method for eigenvalues/eigenvectors)
    vector<vector<double>> eigenvectors(m, vector<double>(m, 0));  // For storing eigenvectors

    // For simplicity, using identity matrix as initial guess for eigenvectors
    for (int i = 0; i < m; i++) {
        eigenvectors[i][i] = 1;
    }

    // Eigenvalues (diagonal values of the covariance matrix)
    vector<double> eigenvalues(m, 0);
    for (int i = 0; i < m; i++) {
        eigenvalues[i] = cov_matrix[i][i];  // Simplified for this demonstration
    }

    // Step 4: Sort eigenvalues and corresponding eigenvectors in descending order
    vector<pair<double, vector<double>>> eigen_pairs;
    for (int i = 0; i < m; i++) {
        eigen_pairs.push_back({eigenvalues[i], eigenvectors[i]});
    }

    sort(eigen_pairs.begin(), eigen_pairs.end(), greater<pair<double, vector<double>>>());

    cout << "Eigenvectors computed..." << "\n";

    // Select the top 'num_components' eigenvectors (principal components)
    vector<vector<double>> top_eigenvectors;
    for (int i = 0; i < num_components; i++) {
        top_eigenvectors.push_back(eigen_pairs[i].second);
    }

    // Step 5: Project the data onto the top eigenvectors
    vector<vector<double>> transformed_data(n, vector<double>(num_components, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < num_components; j++) {
            for (int k = 0; k < m; k++) {
                transformed_data[i][j] += centered[i][k] * top_eigenvectors[j][k];
            }
        }
    }

    cout << "Projection completed..." << "\n";

    // Return the eigenvectors, eigenvalues, and the transformed data
    return make_tuple(top_eigenvectors, eigenvalues, transformed_data);

    cout << "PCA Computation Finished.";
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
        // Optionally, you can store or process the header here
    }

    while (getline(file, line)) { // Read each remaining line from the file
        stringstream ss(line); // Create a stringstream object for splitting the line
        string cell;
        vector<double> row; // To store the split cells

        bool firstCell = true; // Flag to skip the first cell
        while (getline(ss, cell, ',')) { // Split by comma
            if (firstCell) {
                firstCell = false; // Skip the first cell in the row
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

    

    int num_components = 2;  // Number of principal components to keep
    // pca(data, num_components);
    auto [eigenvectors, eigenvalues, pca_result] = pca(data, num_components);

    // Output the PCA result
    cout << "Eigenvectors:\n";
    for (const auto& vec : eigenvectors) {
        for (double val : vec) {
            cout << val << " ";
        }
        cout << endl;
    }

    cout << "\nEigenvalues:\n";
    for (double val : eigenvalues) {
        cout << val << " ";
    }
    cout << endl;

    cout << "\nPCA Transformed Data (First " << num_components << " Components):\n";
    for (const auto& row : pca_result) {
        for (double value : row) {
            cout << value << " ";
        }
        cout << endl;
    }

    return 0;
}