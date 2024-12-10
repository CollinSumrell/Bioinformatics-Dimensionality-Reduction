# PCA Visualizer

**Authors:** Charlie Street, Collin Sumrell
**Affiliation:** University of Oklahoma, Gallogly College of Engineering, School of Computer Science

## Overview

PCA Visualizer is a project aimed at implementing Principal Component Analysis (PCA) from scratch, complemented by an educational visualization tool. The tool provides animations to help users intuitively understand PCA and visualize clustering data.

This project also includes utilities for preprocessing datasets, clustering, and generating animations that showcase the mechanics and results of PCA.

---

## Features

- Custom PCA Implementation: PCA from scratch in Python and C++ for educational purposes, avoiding pre-built library functions like numpy.linalg.eig or scikit-learn.
- Visualization: Animations demonstrating PCA and dimensionality reduction in 2D and 3D.
- Clustering Integration: Support for visualizing clustering data alongside PCA results, with a utility to generate clustering assignments using K-means.
- Data Utilities: Tools to preprocess and convert datasets, particularly for transforming sparse matrices to a more usable CSV format.
Project Motivation

---

## Project Motivation

Bioinformatics and other fields often work with high-dimensional datasets where dimensionality reduction techniques like PCA are essential for analysis and visualization.
This project was developed to:

- Provide a deeper understanding of PCA mechanics through a from-scratch implementation.
- Offer educational animations to make PCA accessible to students and researchers.
- Create a customizable and intuitive tool for exploring PCA results.
Requirements

---

## Requirements

Python: Version 3.8+
- Required Libraries: numpy, scikit-learn (for clustering), pandas
- C++ Compiler: Supporting C++17 or higher
- R (optional): For data preprocessing
Installation

---

## Installation

Clone the repository:
git clone https://github.com/yourusername/pca-visualizer.git  
cd pca-visualizer  
Install dependencies (Python):
pip install -r requirements.txt  
Build the C++ PCA implementation (if using):
cd cpp/  
make  
Usage

--- 

## Running PCA
Preprocess your dataset if necessary (see /utils/preprocessing for scripts).
Run the Python or C++ PCA implementation on your dataset:
python pca.py --input data.csv --output pca_results.csv  
./cpp/pca --input data.csv --output pca_results.csv  
Generating Animations
Use the clustering utility to generate clustering data if needed:
python cluster.py --input data.csv --clusters 5 --output clusters.csv  
Generate an animation from the results:
python visualize.py --pca pca_results.csv --clusters clusters.csv  

MIT License

