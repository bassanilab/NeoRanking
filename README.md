# README (under construction!!!!)
Priorization of neoantigens by machine learning

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
