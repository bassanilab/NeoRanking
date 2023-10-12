# test LR and XGBoost classifiers on NCI_test, TESLA and HiTIDE for mutation and neopep data

source configure.sh

# test logistic regression
cmd="PYTHONPATH=$NEORANKING_CODE python3 Classifier/TestClassifiers.py -c LR_*short*.sav -tr NCI_train -te NCI_test -te TESLA -te HiTIDE -pt neopep"
echo $cmd
eval $cmd

cmd="PYTHONPATH=$NEORANKING_CODE python3 Classifier/TestClassifiers.py -c LR_*long*.sav -tr NCI_train -te NCI_test -te TESLA -te HiTIDE -pt mutation"
echo $cmd
eval $cmd

# test XGBoost
cmd="PYTHONPATH=$NEORANKING_CODE python3 Classifier/TestClassifiers.py -c XGBoost_*short*.xgbm -tr NCI_train -te NCI_test -te TESLA -te HiTIDE -pt neopep"
echo $cmd
eval $cmd

cmd="PYTHONPATH=$NEORANKING_CODE python3 Classifier/TestClassifiers.py -c XGBoost_*long*.xgbm -tr NCI_train -te NCI_test -te TESLA -te HiTIDE -pt mutation"
echo $cmd
eval $cmd
