# Priorization of neoantigens by machine learning

## Introduction
This python code was written to rank neoantigens according to the probability that they are recognized by CD8 T cells. Large data matrices consisting of thousands of neoantigens from 131 cancer patients annotated with several feature scores are used to train machine learning classifiers that rank the neoantigens in a test set. Here we provide the python code and shell scripts that preprocess these data matrices, perform classifier training and testing, and plot the figures of our paper [[1](#Citation)]

### Installation for linux

Install python (here we used python version 3.8.10) with the dependencies outlined in the [requirements.txt](https://github.com/bassanilab/NeoRanking/blob/master/requirements.txt) file:
```
pip install -r requirements.txt
```
Edit the configure.sh file and set the environment variables NEORANKING_RESOURCE for the data directory and NEORANKING_CODE for the code directory. Source the configure.sh file:
```
source configure.sh
```
This will create the data and code directories and various subdirectories if these directories do not yet exist. Download the python code from this github repository (https://github.com/bassanilab/NeoRanking.git) and place it into the $NEORANKING_CODE directory. Download the data matrices from the links indicated here or in [[1](#Citation)] and place the files [Mutation_data_org.txt](https://figshare.com/s/3c27fa3b705a74bdfa10) and [Neopep_data_org.txt](https://figshare.com/s/a000b0990465ab3e9d33) into the $NEORANKING_RESOURCE/data directory, and HLA_allotypes.txt into the $NEORANKING_RESOURCE/hla directory.

If you wish to recreate the plots for Figures 1B, S2A-C, in the paper you need to download the MmpsTestingSet.txt, MmpsTrainingSet.txt, NmersTestingSet.txt, and NmersTrainingSet.txt files from the figshare links provided by Gartner et al. [[2](#Citation)]. These files contain mutations (nmers) and neo-peptides (mmps) together with feature scores and immunogenicity screening annotations used by Gartner et al. If you wish to recreate Figures 3D-F you need to download the file mmc5.xlsx from the Supplemental Data in Wells et al. [[3](#Citation)]

### Running the code

1) Preprocess the original data matrices Mutation_data_org.txt and Neopep_data_org.txt (necessary preprocessing step to be performed once at the start of the analysis). Preprocessing consists of several steps: a) Select the SNV mutations. b) Calculate numerical encoding values for categorical features. c) Impute missing values. d) Transform values of numerical features by quantile normalization. e) Replace categories by encoded numerical values.
    ```
    bash preprocess_data.sh
    ```
2) Training the classifiers. This is only required if training needs to be done on different data or repeated with different parameters. Otherwise classifier models for logistic regression and XGBoost trained on NCI-train [[1](#Citation)] can be obtained from figshare for [neo-peptides](https://figshare.com/s/a000b0990465ab3e9d33) and [mutations](https://figshare.com/s/3c27fa3b705a74bdfa10):
    ```
    bash train_classifier.sh
    ```

3) Testing the classifiers: 
    ```
    bash test_classifier.sh
    ```
    Precalculated [classifier result files](https://figshare.com/s/9fc6c11691273efe995e) used in [[1](#Citation)] can be downloaded from figshare. Place them in the ```classifier_results``` directory and respective subdirectories.

4) Plot figure X: 
    ```
    bash plot_figure_X.sh
    ```
    The plots in Figures 3 and 4, and Suppl. Figure 4 require classifier result files (see above). If you want to reproduce the figures from the paper [[1](#Citation)] based on the results presented there, you can run the scripts as they are. If you prefer to train your own classifiers and plot the figures based on these results, please adapt the corresponding paths and regular expressions in the scripts. Some plots may look slightly different from the ones in the paper due to the random components (especially Shapley values are subject to variations). If you retrained the classifiers, there might be small differences due to the random sampling of non-immunogenic neo-peptides and the stochastic hyperopt parameter optimization.

### Licence
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland

### Contact

For questions regarding code and machine learning methods, please contact Markus Müller (markus.muller@chuv.ch)

For any other questions, please contact Michal Bassani-Sternberg (michal.bassani@chuv.ch)

### Citation

1. Müller M, Huber F, Arnaud M, Kraemer A, Ricart Altimiras E, Michaux J, Taillandier-Coindard M, Chiffelle J, Murgues B, Gehret T, Auger A, Stevenson BJ, Coukos G, Harari A, Bassani-Sternberg M
Machine learning methods and harmonized datasets improve immunogenic neoantigen prediction, in press
2. Gartner JJ, Parkhurst MR, Gros A, Tran E, Jafferji MS, Copeland A, Hanada K-I, Zacharakis N, Lalani A, Krishna S, et al. (2021). A machine learning model for ranking candidate HLA class I neoantigens based on known neoepitopes from multiple human tumor types. Nat. Cancer, 1–12. 10.1038/s43018-021-00197-6
3. Wells DK, van Buuren MM, Dang KK, Hubbard-Lucey VM, Sheehan KCF, Campbell KM, Lamb A, Ward JP, Sidney J, Blazquez AB, et al. (2020). Key Parameters of Tumor Epitope Immunogenicity Revealed Through a Consortium Approach Improve Neoantigen Prediction. Cell 183, 818-834.e13. 10.1016/j.cell.2020.09.015



