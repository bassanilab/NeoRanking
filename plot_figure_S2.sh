source configure.sh

cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeaturePairs.py -fn Figure_S2A -pt neopep -ds NCI -ft pdf -fpair mutant_rank,TAP_score -fih 25 -fiw 25 -tis 13 -las 20 -les 15 -rot 40 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeaturePairs.py -fn Figure_S2B -pt neopep -ds NCI -ft pdf -fpair mutant_rank,mut_netchop_score_ct -fih 25 -fiw 25 -tis 13 -las 20 -les 15 -rot 40 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeaturePairs.py -fn Figure_S2C -pt neopep -ds NCI -ft pdf -fpair mutant_rank,mut_Rank_Stab -fih 25 -fiw 25 -tis 13 -las 20 -les 15 -rot 40 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeatureHistos.py -fn Figure_S2D -pt neopep -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f mutant_rank -fih 5 -fiw 5 -rot 45 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeatureHistos.py -fn Figure_S2E -pt mutation -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f COUNT_MUT_RANK_CI_MIXMHC -fih 5 -fiw 5 -rot 0 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeatureHistos.py -fn Figure_S2F -pt neopep -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f mutant_rank_PRIME -fih 5 -fiw 5 -rot 45 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeatureHistos.py -fn Figure_S2G -pt neopep -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f DAI_MixMHC -fih 5 -fiw 5 -rot 0 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeatureHistos.py -fn Figure_S2H -pt neopep -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f mut_is_binding_pos -fih 5 -fiw 5 -rot 0 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeatureHistos.py -fn Figure_S2I -pt neopep -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f pep_mut_start_10 -fih 5 -fiw 5 -rot 0 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeatureHistos.py -fn Figure_S2J -pt neopep -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f seq_len -fih 5 -fiw 5 -rot 0 -nt -o 8 9 10 11 12"
echo $cmd
eval $cmd
