source configure.sh

CLF=$1  # one of: SVM SVM-lin RF CART ADA LR NNN XGBoost CatBoost TabNet DNN XGBoost
PEPTIDE_TYPE=$2
FEATURE_TYPE=$3
NR_ITER="${4:-200}"
NR_REP="${5:-10}"
NR_NEGATIVE="${6:-100000}"
ALPHA="${7:-0.005}"
NORMALIZER="${8:-q}"
TRAIN_SET="${9:-$(echo $TRAIN_SET)}"
TRAIN_SET_TAG=$(echo $TRAIN_SET | tr '/' '_')
if [ $PEPTIDE_TYPE == 'short' ] 
then
	INPUT_FILE_TAG="${10:-$(echo $INPUT_FILE_TAG_SHORT)}"
else 
	INPUT_FILE_TAG="${10:-$(echo $INPUT_FILE_TAG_LONG)}"
fi
OUT_FILE_TAG="train_${FEATURE_TYPE}_${TRAIN_SET_TAG}_${NR_ITER}_${NR_NEGATIVE}_${ALPHA}_${NORMALIZER}"
OUT_DIR=/home/localadmin/Priorization/Classifiers/

if [ $PEPTIDE_TYPE == 'short' ] && [ $FEATURE_TYPE == 'all' ]
then
	FEATURES=$FEATURE_SHORT_ALL
elif [ $PEPTIDE_TYPE == 'short' ] && [ $FEATURE_TYPE == 'sel' ]
then
	FEATURES=$FEATURE_SHORT_ALL_SEL
elif [ $PEPTIDE_TYPE == 'short' ] && [ $FEATURE_TYPE == 'noMS' ]
then
	FEATURES=$FEATURE_SHORT_NO_IPMSDB
elif [ $PEPTIDE_TYPE == 'short' ] && [ $FEATURE_TYPE == 'noINT' ]
then
	FEATURES=$FEATURE_SHORT_NO_INTOGEN
elif [ $PEPTIDE_TYPE == 'short' ] && [ $FEATURE_TYPE == 'noMSINT' ]
then
	FEATURES=$FEATURE_SHORT_NO_IPMSDB_INTOGEN
elif [ $PEPTIDE_TYPE == 'short' ] && [ $FEATURE_TYPE == 'set1' ]
then
	FEATURES=$FEATURE_SHORT_SET1
elif [ $PEPTIDE_TYPE == 'short' ] && [ $FEATURE_TYPE == 'set2' ]
then
	FEATURES=$FEATURE_SHORT_SET2
elif [ $PEPTIDE_TYPE == 'short' ] && [ $FEATURE_TYPE == 'set3' ]
then
	FEATURES=$FEATURE_SHORT_SET3
elif [ $PEPTIDE_TYPE == 'short' ] && [ $FEATURE_TYPE == 'set4' ]
then
	FEATURES=$FEATURE_SHORT_SET4
elif [ $PEPTIDE_TYPE == 'short' ] && [ $FEATURE_TYPE == 'set5' ]
then
	FEATURES=$FEATURE_SHORT_SET5
elif [ $PEPTIDE_TYPE == 'short' ] && [ $FEATURE_TYPE == 'set6' ]
then
	FEATURES=$FEATURE_SHORT_SET6
elif [ $PEPTIDE_TYPE == 'short' ] && [ $FEATURE_TYPE == 'set7' ]
then
	FEATURES=$FEATURE_SHORT_SET7
elif [ $PEPTIDE_TYPE == 'short' ] && [ $FEATURE_TYPE == 'set8' ]
then
	FEATURES=$FEATURE_SHORT_SET8
elif [ $PEPTIDE_TYPE == 'short' ] && [ $FEATURE_TYPE == 'set9' ]
then
	FEATURES=$FEATURE_SHORT_SET9
elif [ $PEPTIDE_TYPE == 'short' ] && [ $FEATURE_TYPE == 'set10' ]
then
	FEATURES=$FEATURE_SHORT_SET10
elif [ $PEPTIDE_TYPE == 'short' ] && [ $FEATURE_TYPE == 'set11' ]
then
	FEATURES=$FEATURE_SHORT_SET11
elif [ $PEPTIDE_TYPE == 'short' ] && [ $FEATURE_TYPE == 'set12' ]
then
	FEATURES=$FEATURE_SHORT_SET12
elif [ $PEPTIDE_TYPE == 'short' ] && [ $FEATURE_TYPE == 'set13' ]
then
	FEATURES=$FEATURE_SHORT_SET13
elif [ $PEPTIDE_TYPE == 'short' ] && [ $FEATURE_TYPE == 'mut' ]
then
	FEATURES=$FEATURE_SHORT_MUT_SEL
elif [ $PEPTIDE_TYPE == 'long' ] && [ $FEATURE_TYPE == 'all' ]
then
	FEATURES=$FEATURE_LONG_ALL
elif [ $PEPTIDE_TYPE == 'long' ] && [ $FEATURE_TYPE == 'sel' ]
then
	FEATURES=$FEATURE_LONG_ALL_SEL
elif [ $PEPTIDE_TYPE == 'long' ] && [ $FEATURE_TYPE == 'mut' ]
then
	FEATURES=$FEATURE_LONG_MUT_SEL
elif [ $PEPTIDE_TYPE == 'long' ] && [ $FEATURE_TYPE == 'noMS' ]
then
	FEATURES=$FEATURE_LONG_NO_IPMSDB
elif [ $PEPTIDE_TYPE == 'long' ] && [ $FEATURE_TYPE == 'noINT' ]
then
	FEATURES=$FEATURE_LONG_NO_INTOGEN
elif [ $PEPTIDE_TYPE == 'long' ] && [ $FEATURE_TYPE == 'noMSINT' ]
then
	FEATURES=$FEATURE_LONG_NO_IPMSDB_INTOGEN
fi

#for CLF in SVM SVM-lin RF CART ADA LR NNN XGBoost CatBoost TabNet DNN XGBoost
if [ $CLF == 'LR' ] || [ $CLF == 'SVM' ] || [ $CLF == 'SVM-lin' ] || [ $CLF == 'NNN' ] || [ $CLF == 'TabNet'  ] || [ $CLF == 'XGBoost'  ] 
then
	CAT_TO_NUM="-cat float"
elif [ $CLF == 'CatBoost' ]
then
	CAT_TO_NUM="-cat int"
else
	CAT_TO_NUM=""
fi	

if [ $TRAIN_SET == 'HiTIDE' ] || [ $TRAIN_SET == 'TESLA' ]
then
	RESPONSE_TYPES=$RESPONSE_TYPES_TEST
else
	RESPONSE_TYPES=$RESPONSE_TYPES_TRAIN
fi

for i in $(seq 1 $NR_REP) 
do
	cmd="PYTHONPATH=/home/localadmin/Priorization/python python3 Classifier/TrainClassifier.py -c $CLF -s $SCORER -tr $TRAIN_SET -i $INPUT_FILE_TAG -id $OUT_FILE_TAG -f $FEATURES -n $NORMALIZER -v 2 -ni $NR_ITER -cv $NR_CV -a $ALPHA  -mt SNV -rt $RESPONSE_TYPES -im CD8 CD4/CD8 $SHUFFLE $CAT_TO_NUM -pt $PEPTIDE_TYPE -nn $NR_NEGATIVE -eg $EXCL_GENES &"
	echo $cmd
	eval $cmd
	sleep 2s
done


