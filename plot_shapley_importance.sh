source configure.sh
CLF_DIR=$1
CLF_WC1=$2  
CLF_WC2=$3  
DATASET1=$4
DATASET2=$5
ENCODING=$6
PLOT_FILE_PREFIX=$7
PEPTIDE_TYPE=$8
NORMALIZER="${9:-q}"
NR_NEGATIVES="${10:-100000}"
ALPHA="${11:-0.005}"
COLOR_MAP=${12}

if [ "$PEPTIDE_TYPE" == 'short' ] 
then
	FEATURES=$FEATURE_SHORT_ALL
	INPUT_FILE_TAG=$INPUT_FILE_TAG_SHORT
elif [ "$PEPTIDE_TYPE" == 'long' ] 
then
	FEATURES=$FEATURE_LONG_ALL
	INPUT_FILE_TAG=$INPUT_FILE_TAG_LONG
	NR_NEGATIVE=-1
fi

if [ "$ENCODING" == 'float' ]; 
then
	CAT_TO_NUM="-cat float"
elif [ "$ENCODING" == 'int' ];
then
	CAT_TO_NUM="-cat int"
else
	CAT_TO_NUM=""
fi

echo "---> $COLOR_MAP"
if [ -z $COLOR_MAP ] & [ "$COLOR_MAP" != "" ];
then
	COLOR_MAP="-cm \"${COLOR_MAP}\""
fi	

for i in $(seq 1 $NR_REP) 
do
	cmd="PYTHONPATH=/home/localadmin/Priorization/python python3 Scripts/paper/PlotShapleyImportance.py -d $CLF_DIR -c1 $CLF_WC1 -c2 $CLF_WC2 -cat $ENCODING -fp $PLOT_FILE_PREFIX -s $SCORER -ds1 $DATASET1 -ds2 $DATASET2 -i $INPUT_FILE_TAG -f $FEATURES -n $NORMALIZER -v 2 -a $ALPHA  -mt SNV -im CD8 CD4/CD8 $CAT_TO_NUM -pt $PEPTIDE_TYPE -nn $NR_NEGATIVES -a $ALPHA -fd $FEATURE_PLOT_NAMES -las 30 -tis 23 -tts 25 -les 20 $COLOR_MAP &"

	echo $cmd
	eval $cmd
	sleep 2s
done

