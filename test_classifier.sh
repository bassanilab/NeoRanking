source configure.sh

cmd="PYTHONPATH=$NEORANKING_CODE python3 Classifier/TestClassifiers.py -c LR_paper_data*neopep*.sav -tr NCI_train -te NCI_test -te TESLA -te HiTIDE -pt neopep"
echo $cmd
eval $cmd

cmd="PYTHONPATH=$NEORANKING_CODE python3 Classifier/TestClassifiers.py -c LR_paper_data*mutation*.sav -tr NCI_train -te NCI_test -te TESLA -te HiTIDE -pt mutation"
echo $cmd
eval $cmd
