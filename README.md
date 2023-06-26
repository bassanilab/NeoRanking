# Priorization of neoantigens by machine learning

## Introduction
The accurate selection of neoantigens that bind to class I human leukocyte antigen (HLA) and are recognized by autologous T cells is a crucial step in many cancer immunotherapy pipelines.
Here we provide the python code that preprocesses the data matrices, performs classifier training and testing, and plots the figures of our paper [1]

## Installation

Install python with dependencied outlined in the python_dependency.txt file. Download the python code

## Contact

For questions regarding code and machine learning methods, please contact Markus Müller (markus.muller@chuv.ch)
For any other questions, please contact Michal Bassani-Sternberg (michal.bassani@chuv.ch)

## Citation

1. Müller M, Huber F, Arnaud M, Kraemer A, Ricart Altimiras E, Michaux J, Taillandier-Coindard M, Chiffelle J, Murgues B, Gehret T, Auger A, Stevenson BJ, Coukos G, Harari A, Bassani-Sternberg M
Machine learning methods and harmonized datasets improve immunogenic neoantigen prediction, Under revision

1) Download the mutation and neo-peptide data matrices from figshare as described in the paper
2) Edit the paths in GlobalParameters.py file to match your environment
3) Select the row that are used for machine learning (ML). These are the rows that contain single nucleotide variants (SNV)
   and some rows with missing values are removed.
   bash select_ml_data.sh
4) Calculate the categorical encodings for the ML data:
   bash calc_cat_encodings.sh
5) Perform missing value imputation, data normalization and apply categorical feature encoding. Store the normalized data files:
   bash normalize_data.sh
6) Plot Figure 1B-I:
   bash plot_figure_1.sh
7) Plot Figure 2:
   bash plot_figure_2.sh
...
)
