source configure.sh

# select rows used for ML. These are all somatic SNV mutations.
#PYTHONPATH=$NEORANKING_CODE python3 DataWrangling/SelectMLData.py -pt neopep
#PYTHONPATH=$NEORANKING_CODE python3 DataWrangling/SelectMLData.py -pt mutation

# calculate numerical encoding for categorical values based on NCI_train
#PYTHONPATH=$NEORANKING_CODE python3 DataWrangling/CalcCatEncodings.py -pt neopep -ds NCI_train
#PYTHONPATH=$NEORANKING_CODE python3 DataWrangling/CalcCatEncodings.py -pt mutation -ds NCI_train

# normalize (missing value imputation, data normalization, and replacing categories by their encoded numerical values
PYTHONPATH=$NEORANKING_CODE python3 DataWrangling/NormalizeData.py -pt neopep -ds NCI_train -o ml
PYTHONPATH=$NEORANKING_CODE python3 DataWrangling/NormalizeData.py -pt mutation -ds NCI_train -o ml
PYTHONPATH=$NEORANKING_CODE python3 DataWrangling/NormalizeData.py -pt neopep -ds NCI_train -o plot
PYTHONPATH=$NEORANKING_CODE python3 DataWrangling/NormalizeData.py -pt mutation -ds NCI_train -o plot
