PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PrintDatasetPatients.py -pt long
PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PrintDatasetPatients.py -pt short

PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PrintDatasetStats.py -i rt_stab_chop_tap_mbp_tcr -pt short -mt SNV -rt CD8 CD4/CD8 negative not_tested -im CD8 CD4/CD8
PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PrintDatasetStats.py -i rt_netmhc_stab_chop_tap_mbp_tcr -pt long -mt SNV -rt CD8 CD4/CD8 negative not_tested -im CD8 CD4/CD8
