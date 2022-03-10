for neodisc_file in $(ls /home/localadmin/Priorization/results/*rt_netmhc_stab.txt)
do
	new_file=(`echo $neodisc_file | sed 's/rt/long_rt/g'`)

	if [ ! -f "$new_file" ]; then
		eval "mv $neodisc_file $new_file"
	else
		echo "File $new_file exists."
	fi
done

