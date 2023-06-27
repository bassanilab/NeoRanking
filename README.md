# Priorization of neoantigens by machine learning

## Introduction
This python code was written to rank neoantigens according to their probability that they are recognized by CD8 T cells. Large data matrices consisting of thousands neoantigens from 131 cancer patients annotated with several feature scores are used to train machine learning classifiers that rank the neoantigens in a test set. Here we provide the python code and shell scripts that preprocess these data matrices, perform classifier training and testing, and plot the figures of our paper [[1](#Citation)]

### Installation for linux

Install python with the dependencies outlined in the [requirements.txt](https://github.com/bassanilab/NeoRanking/blob/master/requirements.txt) file:
```
pip install -r requirements.txt
```
Edit the configure.sh file and set the environment variables NEORANKING_RESOURCE for the data directory and NEORANKING_CODE for the code directory. Source the configure.sh file:
```
source configure.sh
```
This will create the data and code directories and various subdirectories, if these directories do not yet exist. Download the python code from this github repository (https://github.com/bassanilab/NeoRanking.git) and place it into the $NEORANKING_CODE directory. Download the data matrices from the links indicated in [[1](#Citation)] and place the files Mutation_data_org.txt and Neopep_data_org.txt into the $NEORANKING_RESOURCE/data directory, and HLA_allotypes.txt into the $NEORANKING_RESOURCE/hla directory.

### Running the code

1) Preprocess the original data matrices Mutation_data_org.txt and Neopep_data_org.txt (necessary preprocessing step to be performed once at the start of the analysis). Preprocessing consists of several steps: 1) Select the SNV mutations. 2) Calculalate numerical encoding values for categorical features. 3) Impute missing values. 4) Transform values of numerical features by quantile normalization. 5) Replace cetagoricies by encoded numerical vealues.
```
bash preprocess_data.sh
```
2) Training the classifiers. This is only required if training needs to be done on different data or repeated with different parameters. Otherwise pretrained classifiers can be obtained from the links in [[1](#Citation)]: 
```
bash train_classifier.sh
```
3) Testing the classifiers: 
```
bash test_classifier.sh
```
4) Plot figure X: 
```
bash plot_figure_X.sh
```
### Licence

Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland

### Contact

For questions regarding code and machine learning methods, please contact Markus Müller (markus.muller@chuv.ch)

For any other questions, please contact Michal Bassani-Sternberg (michal.bassani@chuv.ch)

### Citation

1. Müller M, Huber F, Arnaud M, Kraemer A, Ricart Altimiras E, Michaux J, Taillandier-Coindard M, Chiffelle J, Murgues B, Gehret T, Auger A, Stevenson BJ, Coukos G, Harari A, Bassani-Sternberg M
Machine learning methods and harmonized datasets improve immunogenic neoantigen prediction, Under revision

