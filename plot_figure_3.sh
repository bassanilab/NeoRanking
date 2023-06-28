source configure.sh

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotTESLAMLVotingComp.py -fn Figure_3D -ft pdf -pt FR_plot"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotTESLAMLVotingComp.py -fn Figure_3E -ft pdf -pt TTIF_plot"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotTESLAMLVotingComp.py -fn Figure_3F -ft pdf -pt AUPRC_plot"
echo $CMD
eval $CMD
