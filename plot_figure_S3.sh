source configure.sh

cmd="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 Scripts/paper/PlotFeaturePairs.py -fn Suppl_Figure_3A -pt neopep -ds NCI -ft pdf -fpair mutant_rank,mut_netchop_score_ct -fih 25 -fiw 25 -tis 13 -las 20 -les 15 -rot 40 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 Scripts/paper/PlotFeaturePairs.py -fn Suppl_Figure_3B -pt neopep -ds NCI -ft pdf -fpair mutant_rank,TAP_score -fih 25 -fiw 25 -tis 13 -las 20 -les 15 -rot 40 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 Scripts/paper/PlotFeaturePairs.py -fn Suppl_Figure_3C -pt neopep -ds NCI -ft pdf -fpair mutant_rank,mut_Rank_Stab -fih 25 -fiw 25 -tis 13 -las 20 -les 15 -rot 40 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 Scripts/paper/PlotFeatureHistos.py -fn Suppl_Figure_3D -pt neopep -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f mutant_rank -fih 5 -fiw 5 -rot 45 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 Scripts/paper/PlotFeatureHistos.py -fn Suppl_Figure_3E -pt mutation -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f COUNT_MUT_RANK_CI_MIXMHC -fih 5 -fiw 5 -rot 0 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 Scripts/paper/PlotFeatureHistos.py -fn Suppl_Figure_3F -pt neopep -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f mutant_rank_PRIME -fih 5 -fiw 5 -rot 45 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 Scripts/paper/PlotFeatureHistos.py -fn Suppl_Figure_3G -pt neopep -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f DAI_MixMHC -fih 5 -fiw 5 -rot 0 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 Scripts/paper/PlotFeatureHistos.py -fn Suppl_Figure_3H -pt neopep -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f mut_is_binding_pos -fih 5 -fiw 5 -rot 0 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 Scripts/paper/PlotFeatureHistos.py -fn Suppl_Figure_3I -pt neopep -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f pep_mut_start_10 -fih 5 -fiw 5 -rot 0 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 Scripts/paper/PlotFeatureHistos.py -fn Suppl_Figure_3J -pt neopep -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f seq_len -fih 5 -fiw 5 -rot 0 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 Scripts/paper/PlotFeatureHistos.py -fn Suppl_Figure_3H -pt mutation -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f rnaseq_TPM -fih 5 -fiw 5 -rot 45 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 Scripts/paper/PlotFeatureHistos.py -fn Suppl_Figure_3I -pt mutation -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f rnaseq_alt_support -fih 5 -fiw 5 -rot 45 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 Scripts/paper/PlotFeatureHistos.py -fn Suppl_Figure_3J -pt mutation -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f bestWTPeptideCount_I -fih 5 -fiw 5 -rot 45 -nt"
echo $cmd
eval $cmd
