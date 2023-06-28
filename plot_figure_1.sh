source configure.sh

if [ ! -f "$NEORANKING_RESOURCE/data/Compare_in-house_Gartner_Counts.txt" ]; then
  CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/CompareInHouseGartnerCounts.py"
  echo $CMD
  eval $CMD
fi

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotInHouseGartnerMutationCounts.py -fn Figure_1B -pt SNV_count -ft pdf -las 50 -tis 25 -les 15 -ps 300  -fiw 8 -fih 8 -cmp viridis"
echo $CMD
eval $CMD

if [[ ! -f "$NEORANKING_RESOURCE/data/Patient_statistics_neopep.txt" || ! -f "$NEORANKING_RESOURCE/data/Patient_statistics_mutation.txt" ]]; then
  CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PrintDatasetStats.py"
  echo $CMD
  eval $CMD
fi

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotDatasetStats.py -fn Figure_1C -pt mutation -s \"Mut-seq count\" -ft pdf -rot 0 -las 40 -tis 30 -les 20 -fiw 10 -fih 6"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotDatasetStats.py -fn Figure_1D -pt neopep -s \"Neo-pep count\" -ft pdf -rot 0 -las 40 -tis 30 -les 20 -fiw 10 -fih 6"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotDatasetStats.py -fn Figure_1E -pt neopep -s \"MixMHCpred rank\" -ft pdf -rot 0 -las 40 -tis 30 -les 20 -fiw 10 -fih 6"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotDatasetStats.py -fn Figure_1F -pt mutation -s \"RNAseq TPM\" -ft pdf -rot 0 -las 40 -tis 30 -les 20 -fiw 10 -fih 6"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotDatasetStats.py -fn Figure_1G -pt mutation -s \"RNAseq coverage\" -ft pdf -rot 0 -las 40 -tis 30 -les 20 -fiw 10 -fih 6"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotDatasetStats.py -fn Figure_1H -pt neopep -s \"Neo-pep_imm count per mut-seq\" -ft pdf -rot 0 -las 40 -tis 30 -les 20 -fiw 10 -fih 6"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotCD8MutationsCnt.py -fn Figure_1I -ft pdf -fih 20 -fiw 30 -las 25 -tis 20 -les 20 -ps 100"
echo $CMD
eval $CMD

