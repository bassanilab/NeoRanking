source configure.sh

if [ ! -f "$NEORANKING_RESOURCE/data/Compare_in-house_Gartner_Counts.txt" ]; then
  CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/CompareInHouseGartnerCounts.py"
  echo $CMD
  eval $CMD
fi

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotInHouseGartnerMutationCounts.py -fn Figure_S2A -pt SNV_imm_count -ft pdf -las 40 -tis 30 -les 20 -fiw 10 -fih 10"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotInHouseGartnerMutationCounts.py -fn Figure_S2B -pt InDel_count -ft pdf  -las 40 -tis 30 -les 20 -fiw 10 -fih 10"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotInHouseGartnerMutationCounts.py -fn Figure_S2C -pt FS_count -ft pdf  -las 40 -tis 30 -les 20 -fiw 10 -fih 10"
echo $CMD
eval $CMD

