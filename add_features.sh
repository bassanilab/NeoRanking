source configure.sh

in_tag=$1
feature_tag=$2
out_tag=$3
peptide_type=$4
for neodisc_file in $(eval "ls /home/localadmin/Priorization/results/*${peptide_type}_${in_tag}.txt")
do
	out_file=(`echo $neodisc_file | sed "s/\.txt/_$out_tag\.txt/g"`)
	fn=${neodisc_file##*/}
	p=(`echo $fn | sed 's/_.*//g'`)

	if [ "$p" != '13LN' ]; then
		continue
	fi

	if [ ! -f "$out_file" ]; then
		if [ "$feature_tag" == "pred" ]; then
			if [ "$peptide_type" == "long" ]; then
				classifier=$CLF_LONG
				features=$FEATURE_LONG_MUT_SEL
			else
				classifier=$CLF_SHORT
				features=$FEATURE_SHORT_BIND_SEL
			fi
			cmd="PYTHONPATH=/home/localadmin/Priorization/python python3 Features/AddFeatures.py -p ${p} --input_file_tag $in_tag --$feature_tag --peptide_type $peptide_type \
				-mt SNV -clf $classifier -cf $features -cn $NORMALIZER $CAT_TO_NUM &"
		else
			cmd="PYTHONPATH=/home/localadmin/Priorization/python python3 Features/AddFeatures.py -p ${p} --input_file_tag $in_tag --$feature_tag --peptide_type $peptide_type &"
		fi
		echo $cmd
		eval $cmd
	else
		echo "$feature_tag already processed for patient $p"
	fi
done

