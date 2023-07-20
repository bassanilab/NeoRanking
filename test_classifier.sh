source configure.sh

# test logistic regression
cmd="PYTHONPATH=$NEORANKING_CODE python3 Classifier/TestClassifiers.py -c LR_paper_data*neopep*.sav -tr NCI_train -te NCI_test -te TESLA -te HiTIDE -pt neopep"
echo $cmd
eval $cmd

cmd="PYTHONPATH=$NEORANKING_CODE python3 Classifier/TestClassifiers.py -c LR_paper_data*mutation*.sav -tr NCI_train -te NCI_test -te TESLA -te HiTIDE -pt mutation"
echo $cmd
eval $cmd

# test XGBoost
cmd="PYTHONPATH=$NEORANKING_CODE python3 Classifier/TestClassifiers.py -c XGBoost_paper_data*neopep*.xgbm -tr NCI_train -te NCI_test -te TESLA -te HiTIDE -pt neopep"
echo $cmd
eval $cmd

cmd="PYTHONPATH=$NEORANKING_CODE python3 Classifier/TestClassifiers.py -c XGBoost_paper_data*mutation*.xgbm -tr NCI_train -te NCI_test -te TESLA -te HiTIDE -pt mutation"
echo $cmd
eval $cmd
