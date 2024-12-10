# PCA Visualizer

**Authors:** Charlie Street, Collin Sumrell
**Affiliation:** University of Oklahoma, Gallogly College of Engineering, School of Computer Science

## Overview

PCA Visualizer is a project aimed at implementing Principal Component Analysis (PCA) from scratch, complemented by an educational visualization tool. The tool provides animations to help users intuitively understand PCA and visualize clustering data.

This project also includes utilities for preprocessing datasets, clustering, and generating animations that showcase the mechanics and results of PCA.

---

## Features

- Custom PCA Implementation: PCA from scratch in Python and C++ for educational purposes, avoiding pre-built library functions like `numpy.linalg.eig` or scikit-learn.
- Visualization: Animations demonstrating PCA and dimensionality reduction in 2D and 3D.
- Clustering Integration: Support for visualizing clustering data alongside PCA results, with a utility to generate clustering assignments using K-means.
- Data Utilities: Tools to preprocess and convert datasets, particularly for transforming sparse matrices to a more usable CSV format.

---

## Project Motivation

Bioinformatics and other fields often work with high-dimensional datasets where dimensionality reduction techniques like PCA are essential for analysis and visualization.
This project was developed to:

- Provide a deeper understanding of PCA mechanics through a from-scratch implementation.
- Offer educational animations to make PCA accessible to students and researchers.
- Create a customizable and intuitive tool for exploring PCA results.

---

## Requirements

**Note:** This has only been tested on macOS 14.4+. If you don't have a Mac available to you, we've pre-rendered and pre-run all animations and PCA on the datasets. You can find where to locate those under the Layout header.

- Python: Version 3.8+
  - Required Libraries: numpy, scikit-learn (for clustering), pandas
- C++ Compiler: Supporting C++17 or higher
- R (optional): For data preprocessing
- Homebrew

---

## Installation

Clone the repository:

```bash
git clone https://github.com/CollinSumrell/Bioinformatics-Dimensionality-Reduction.git
cd pca-visualizer
```

(alternatively, you can just download the ZIP)

Install dependencies:

```bash
pip install -r requirements.txt
brew install eigen
```

Build the C++ PCA implementation (if using):

```bash
make
```

---

## Running PCA

Preprocess your dataset if necessary (see /utils/preprocessing for scripts).
Run the Python or C++ PCA implementation on your dataset:

```bash
make # for C++
```

For the Python implementation, just navigate to that file and run the script.

You'll need to specify the path to your dataset in the C++ and Python files, respectively.

## Generating Animations

First, follow the instructions [here](https://docs.manim.community/en/stable/installation.html) to install Manim.

Use the clustering utility to generate clustering data if needed:

```bash
manim -qh src/animations/<filename>.py
```

`h` in the above terminal command can be replaced with `l`, `m`, `h`, or `k` for 480p, 720p, 1080p, or 4K resolution, respectively.

## Layout

- Source files for the animations can be found under `/src/animations`
  - Pre-rendered versions are under `/media/videos`
- Source files for the PCA implementations are under `/src/cpp_pca` and `/src/python_pca/`
- Python and R utilities to help assemble the datasets, resize them, and merge CSVs can be found under `/src/utils`
- The datasets that we ran PCA on are under `/datasets/csv`
- The results of our C++ PCA implementation can be found under `/results/`
  - A Python script that calculates statistics on the PCA output can also be found here. It was used for troubleshooting during development, primarily around how we were standardizing/normalizing the data/results.
