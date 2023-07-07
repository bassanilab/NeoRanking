source configure.sh

cmd="PYTHONPATH=$NEORANKING_CODE python3 Classifier/TestClassifiers.py -c LR*.sav -tr NCI_train -te NCI_test -te TESLA -te HiTIDE -pt neopep"
echo $cmd
eval $cmd

