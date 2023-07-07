source configure.sh

cmd="PYTHONPATH=$NEORANKING_CODE python3 Classifier/TrainClassifier.py -c LR -tr NCI_train -pt neopep -tag paper_data"
echo $cmd
eval $cmd
