# Priorization of neoantigens by machine learning

## Introduction
This python code was written to prioritize neoantigens for cancer immunotherapy. Large data matrices consisting of thousands neoantigens from 131 cancer patients annotated with several feature scores are used to train machine learning classifiers that rank the neoantigens in a test set. Here we provide the python code and shell scripts that preprocess these data matrices, perform classifier training and testing, and plot the figures of our paper [[1](#Citation)]

### Installation

Install python with the dependencies outlined in the [requirements.txt](https://github.com/bassanilab/NeoRanking/blob/master/requirements.txt) file. Download the python code from this github repository (https://github.com/bassanilab/NeoRanking.git). Adapt the paths indicated in Utils/GlobalParameters.py file to your environment. Download the data matrices from the links indicated in [[1](#Citation)] and place the files Mutation_data_org.txt, Neopep_data_org.txt, and HLA_allotypes.txt into the paths indicated in Utils/GlobalParameters.py.

### Running the code

1) Select the rows in Mutation_data_org.txt and Neopep_data_org.txt that are used for machine learning:
'''
bash select_ml_data.sh
'''
2) Calculate categorical encodings: 
'''
bash calc_cat_encodings.sh
'''

### Licence

Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland

### Contact

For questions regarding code and machine learning methods, please contact Markus Müller (markus.muller@chuv.ch)

For any other questions, please contact Michal Bassani-Sternberg (michal.bassani@chuv.ch)

### Citation

1. Müller M, Huber F, Arnaud M, Kraemer A, Ricart Altimiras E, Michaux J, Taillandier-Coindard M, Chiffelle J, Murgues B, Gehret T, Auger A, Stevenson BJ, Coukos G, Harari A, Bassani-Sternberg M
Machine learning methods and harmonized datasets improve immunogenic neoantigen prediction, Under revision

