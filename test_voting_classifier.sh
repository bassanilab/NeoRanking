source configure.sh

cmd="PYTHONPATH=$NEORANKING_CODE python3 Classifier/TestVotingClassifier.py -c1 LR_paper*_neopep_*.sav -c2 XGBoost_paper*_neopep_*.xgbm -w 0.9 -tr NCI_train -te NCI_test -te TESLA -te HiTIDE -pt neopep"
echo $cmd
eval $cmd

