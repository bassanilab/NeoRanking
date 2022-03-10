for neodisc_file in $(ls /home/localadmin/Priorization/results/*rt_netmhc.txt)
do
	netstab_file=(`echo $neodisc_file | sed 's/rt_netmhc/rt_netmhc_stab/g'`)
	fn=${neodisc_file##*/}
	p=(`echo $fn | sed 's/_.*//g'`)

	if [ ! -f "$netatb_file" ]; then
		cmd="PYTHONPATH=/home/localadmin/Priorization/python python3 Features/AddFeatures.py -p ${p} --input_file_tag "rt_netmhc" --netmhcstab &"
		echo $cmd
		eval $cmd
	else
		echo "Netmhc already processed for patient $p"
	fi
done

