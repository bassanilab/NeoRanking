# test LR and XGBoost voting classifier on NCI_test, TESLA and HiTIDE

source configure.sh

cmd="PYTHONPATH=$NEORANKING_CODE python3 Classifier/TestVotingClassifier.py -c1 LR_paper*_neopep_*.sav -c2 XGBoost_paper*_neopep_*.xgbm -w 0.9 -tr NCI_train -te NCI_test -te TESLA -te HiTIDE -pt neopep"
echo $cmd
eval $cmd

cmd="PYTHONPATH=$NEORANKING_CODE python3 Classifier/TestVotingClassifier.py -c1 LR_paper*_mutation_*.sav -c2 XGBoost_paper*_mutation_*.xgbm -w 0.5 -tr NCI_train -te NCI_test -te TESLA -te HiTIDE -pt mutation"
echo $cmd
eval $cmd

