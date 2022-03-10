in_tag=$1
feature_tag=$2
out_tag=$3
peptide_type=$4
for neodisc_file in $(eval "ls /home/localadmin/Priorization/results/*$in_tag.txt")
do
	out_file=(`echo $neodisc_file | sed "s/\.txt/_$out_tag\.txt/g"`)
	fn=${neodisc_file##*/}
	p=(`echo $fn | sed 's/_.*//g'`)

	if [ ! -f "$out_file" ]; then
		cmd="PYTHONPATH=/home/localadmin/Priorization/python python3 Features/AddFeatures.py -p ${p} --input_file_tag $in_tag --$feature_tag --peptide_type $peptide_type &"
		echo $cmd
		eval $cmd
	else
		echo "$feature_tag already processed for patient $p"
	fi
done

