for neodisc_file in $(ls /home/localadmin/Priorization/results/*rt.txt)
do
	netmhc_file=(`echo $neodisc_file | sed 's/rt/rt_netmhc/g'`)
	fn=${neodisc_file##*/}
	p=(`echo $fn | sed 's/_.*//g'`)

	if [ ! -f "$netmhc_file" ]; then
		cmd="PYTHONPATH=/home/localadmin/Priorization/python python3 Features/AddFeatures.py -p ${p} --input_file_tag "rt" --netmhcpan &"
		echo $cmd
		eval $cmd
	else
		echo "Netmhc already processed for patient $p"
	fi
done

