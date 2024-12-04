#include <iostream>
#include <Eigen/Dense>
#include <vector>

using namespace Eigen;
using namespace std;

// Function to perform PCA
void PCA(const MatrixXd& data) {
    // Step 1: Center the data by subtracting the mean of each column
    MatrixXd centered = data.rowwise() - data.colwise().mean();

    // Step 2: Compute the covariance matrix
    MatrixXd covariance = (centered.transpose() * centered) / double(centered.rows() - 1);

    // Step 3: Compute the eigenvalues and eigenvectors
    SelfAdjointEigenSolver<MatrixXd> solver(covariance);
    if (solver.info() != Success) {
        cerr << "Eigenvalue decomposition failed!" << endl;
        return;
    }

    // Step 4: Eigenvalues and eigenvectors
    VectorXd eigenvalues = solver.eigenvalues().real();
    MatrixXd eigenvectors = solver.eigenvectors().real();

    // Step 5: Sort the eigenvalues and corresponding eigenvectors in descending order
    vector<pair<double, VectorXd>> eig_pairs;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        eig_pairs.push_back(make_pair(eigenvalues[i], eigenvectors.col(i)));
    }
    sort(eig_pairs.rbegin(), eig_pairs.rend());

    // Step 6: Project data onto the new space
    int k = 2;  // Number of principal components to keep (adjustable)
    MatrixXd projection(data.rows(), k);
    for (int i = 0; i < k; ++i) {
        projection.col(i) = centered * eig_pairs[i].second;
    }

    // Output the projected data
    cout << "Projected Data (first two principal components):\n" << projection << endl;
}

int main() {
    // Example data: 5 samples, 3 features (rows: samples, columns: features)
    MatrixXd data(5, 3);
    data << 2.5, 2.4, 3.5,
            0.5, 0.7, 1.0,
            2.2, 2.9, 3.1,
            1.9, 2.2, 2.7,
            3.1, 3.0, 4.2;

    cout << "Original Data:\n" << data << endl;

    // Perform PCA
    PCA(data);

    return 0;
}