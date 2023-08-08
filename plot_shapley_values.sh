source configure.sh
CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotShapleyValues.py -c LR_paper_data*.sav -tr NCI_train -p 4350 -pept KTYQGSYGFRR -pt neopep -ft pdf -fn Figure_3H -las 30 -tis 23 -tts 25 -les 20"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotShapleyValues.py -c LR_paper_data*.sav -tr NCI_train -p 4324 -pept DRNIFRHSVV -pt neopep -ft pdf -fn Figure_3I -las 30 -tis 23 -tts 25 -les 20"
echo $CMD
eval $CMD
