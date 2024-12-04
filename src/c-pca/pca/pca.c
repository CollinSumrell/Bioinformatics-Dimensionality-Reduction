#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ROWS 4
#define COLS 2

// Function to calculate the mean of each column
void calculate_mean(double data[ROWS][COLS], double mean[COLS]) {
    for (int j = 0; j < COLS; j++) {
        mean[j] = 0;
        for (int i = 0; i < ROWS; i++) {
            mean[j] += data[i][j];
        }
        mean[j] /= ROWS;
    }
}

// Function to center the data (subtract mean)
void center_data(double data[ROWS][COLS], double mean[COLS]) {
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            data[i][j] -= mean[j];
        }
    }
}

// Function to calculate the covariance matrix
void calculate_covariance(double data[ROWS][COLS], double cov[COLS][COLS]) {
    for (int i = 0; i < COLS; i++) {
        for (int j = 0; j < COLS; j++) {
            cov[i][j] = 0;
            for (int k = 0; k < ROWS; k++) {
                cov[i][j] += data[k][i] * data[k][j];
            }
            cov[i][j] /= (ROWS - 1);
        }
    }
}

// Function to perform eigen decomposition (simplified for 2x2 matrix)
void eigen_decomposition(double cov[COLS][COLS], double eigenvalues[COLS], double eigenvectors[COLS][COLS]) {
    double a = cov[0][0];
    double b = cov[0][1];
    double c = cov[1][0];
    double d = cov[1][1];

    // Calculate eigenvalues
    double trace = a + d;
    double determinant = a * d - b * c;
    double discriminant = sqrt(trace * trace - 4 * determinant);

    eigenvalues[0] = (trace + discriminant) / 2;
    eigenvalues[1] = (trace - discriminant) / 2;

    // Calculate eigenvectors (normalized)
    if (b != 0) {
        eigenvectors[0][0] = eigenvalues[0] - d;
        eigenvectors[1][0] = b;
        eigenvectors[0][1] = eigenvalues[1] - d;
        eigenvectors[1][1] = b;
    } else {
        eigenvectors[0][0] = 1;
        eigenvectors[1][0] = 0;
        eigenvectors[0][1] = 0;
        eigenvectors[1][1] = 1;
    }

    // Normalize eigenvectors
    for (int i = 0; i < COLS; i++) {
        double norm = sqrt(eigenvectors[0][i] * eigenvectors[0][i] + eigenvectors[1][i] * eigenvectors[1][i]);
        eigenvectors[0][i] /= norm;
        eigenvectors[1][i] /= norm;
    }
}

int main() {
    // Example dataset (4 samples, 2 features)
    double data[ROWS][COLS] = {
        {2.5, 2.4},
        {0.5, 0.7},
        {2.2, 2.9},
        {1.9, 2.2}
    };

    // Step 1: Calculate the mean of each column
    double mean[COLS];
    calculate_mean(data, mean);

    // Step 2: Center the data
    center_data(data, mean);

    // Step 3: Calculate the covariance matrix
    double cov[COLS][COLS];
    calculate_covariance(data, cov);

    // Step 4: Perform eigen decomposition
    double eigenvalues[COLS];
    double eigenvectors[COLS][COLS];
    eigen_decomposition(cov, eigenvalues, eigenvectors);

    // Output results
    printf("Eigenvalues:\n");
    for (int i = 0; i < COLS; i++) {
        printf("%lf\n", eigenvalues[i]);
    }

    printf("\nEigenvectors:\n");
    for (int i = 0; i < COLS; i++) {
        printf("[%lf, %lf]\n", eigenvectors[0][i], eigenvectors[1][i]);
    }

    return 0;
}
