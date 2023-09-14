source configure.sh

cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeaturePairs.py -fn Figure_2A -pt neopep -ds NCI -ft pdf -fpair mutant_rank,mutant_rank_netMHCpan -fih 25 -fiw 25 -tis 13 -las 20 -les 15 -rot 40 -hl 0.5 -vl 0.5 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeatureHistos.py -fn Figure_2B -pt neopep -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f mutant_other_significant_alleles -fih 5 -fiw 5 -rot 0 -log -nt -lbls 1 2 3 4 5 6"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeaturePairs.py -fn Figure_2C -pt neopep -ds NCI -ft pdf -fpair mut_is_binding_pos,DAI_MixMHC -fih 7 -fiw 5 -tis 13 -las 20 -les 15 -rot 0 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeaturePairs.py -fn Figure_2D -pt mutation -ds NCI -ft pdf -fpair TCGA_Cancer_expression,rnaseq_TPM -fih 7 -fiw 7 -tis 13 -las 20 -les 15 -rot 40 -reg -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeaturePairs.py -fn Figure_2E -pt mutation -ds NCI -ft pdf -fpair GTEx_all_tissues_expression_mean,rnaseq_TPM -fih 7 -fiw 7 -tis 13 -las 20 -les 15 -rot 40 -reg -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeaturePairs.py -fn Figure_2F -pt mutation -ds NCI -ft pdf -fpair rnaseq_TPM,bestWTPeptideCount_I -fih 7 -fiw 7 -tis 13 -las 20 -les 15 -rot 40 -reg -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeaturePairs.py -fn Figure_2G -pt neopep -ds NCI -ft pdf -fpair bestWTPeptideCount_I,bestWTMatchScore_I -fih 7 -fiw 7 -tis 13 -las 20 -les 15 -rot 40 -reg -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeatureHistos.py -fn Figure_2H -pt neopep -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f bestWTMatchType_I -fih 5 -fiw 5 -rot 30 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeaturePairs.py -fn Figure_2I -pt neopep -ds NCI -ft pdf -fpair mut_is_binding_pos,bestWTMatchType_I -fih 5 -fiw 5 -tis 13 -las 20 -les 15 -rot 40 -o NONE PARTIAL PARTIAL_MUT COVER INCLUDED EXACT -nt"
echo $cmd
eval $cmd
