CMD="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 /home/localadmin/Priorization/paper_code/Scripts/paper/PlotCD8MutationsCnt.py -fn Figure_1B -ft pdf -fih 20 -fiw 30 -las 25 -tis 20 -les 20 -ps 100"
echo $CMD
eval $CMD

CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotNeoDiscGartnerMutationCounts.py -f /home/localadmin/Priorization/Plots/DatasetStats/Compare_NeoDisc_Gartner_Counts.txt -fp Figure_NeoDisc_Gartner_Counts -las 50 -tis 25 -les 15 -ps 300 -dpi 600 -fiw 8 -fih 8 -cmp viridis"
#echo $CMD
#eval $CMD

CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotDatasetStats.py -d /home/localadmin/Priorization/Plots/DatasetStats -fp Figure_Dataset_Stats -pt short -rot 0 -las 40 -tis 30 -les 20 -dpi 600 -fiw 10 -fih 6"
#echo $CMD
#eval $CMD

CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotDatasetStats.py -d /home/localadmin/Priorization/Plots/DatasetStats -fp Figure_Dataset_Stats -pt long -rot 30 -las 40 -tis 30 -les 20 -dpi 600 -fiw 10 -fih 6"
#echo $CMD
#eval $CMD

