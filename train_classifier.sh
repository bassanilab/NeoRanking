# train LR and XGBoost classifiers on NCI_train for mutation and neopep data

source configure.sh

# train logistic regression
#cmd="PYTHONPATH=$NEORANKING_CODE python3 Classifier/TrainClassifier.py -c LR -tr NCI_train -pt neopep -tag paper_data"
#echo $cmd
#eval $cmd

cmd="PYTHONPATH=$NEORANKING_CODE python3 Classifier/TrainClassifier.py -c LR -tr NCI_train -pt mutation -tag paper_data"
echo $cmd
eval $cmd

# train XGBoost
#cmd="PYTHONPATH=$NEORANKING_CODE python3 Classifier/TrainClassifier.py -c XGBoost -tr NCI_train -pt neopep -tag paper_data"
#echo $cmd
#eval $cmd

cmd="PYTHONPATH=$NEORANKING_CODE python3 Classifier/TrainClassifier.py -c XGBoost -tr NCI_train -pt mutation -tag paper_data"
echo $cmd
eval $cmd

