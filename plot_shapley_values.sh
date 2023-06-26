source configure.sh
CLF_BIN=$1  
ENCODING=$2
PNG_PREFIX=$3
PEPTIDE_TYPE=$4
NORMALIZER="${5:-q}"
NR_NEGATIVES="${6:-100000}"
ALPHA="${7:-0.005}"

if [ $PEPTIDE_TYPE == 'short' ] 
then
	FEATURES=$FEATURE_SHORT_ALL
	INPUT_FILE_TAG=$INPUT_FILE_TAG_SHORT
elif [ $PEPTIDE_TYPE == 'long' ] 
then
	FEATURES=$FEATURE_LONG_ALL
	INPUT_FILE_TAG=$INPUT_FILE_TAG_LONG
	NR_NEGATIVE=-1
fi

if [ $ENCODING == 'float' ] 
then
	CAT_TO_NUM="-cat float"
elif [ $ENCODING == 'int' ]
then
	CAT_TO_NUM="-cat int"
else
	CAT_TO_NUM=""
fi	

for i in $(seq 1 $NR_REP) 
do
	cmd="PYTHONPATH=/home/localadmin/Priorization/python python3 Scripts/paper/PlotShapleyValues.py -c $CLF_BIN -png $PNG_PREFIX -s $SCORER -tr $TRAIN_SET -te $TEST_SET -i $INPUT_FILE_TAG -f $FEATURES -n $NORMALIZER -v 2 -a $ALPHA  -mt SNV -im CD8 CD4/CD8 $CAT_TO_NUM -pt $PEPTIDE_TYPE -nn $NR_NEGATIVES -a $ALPHA -fd $FEATURE_PLOT_NAMES -las 30 -tis 20 -tts 20 &"

	echo $cmd
	eval $cmd
	sleep 2s
done


