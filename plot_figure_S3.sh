source configure.sh

cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeatureHistos.py -fn Figure_S3A -pt mutation -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f rnaseq_TPM -fih 5 -fiw 5 -rot 45 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeatureHistos.py -fn Figure_S3B -pt mutation -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f rnaseq_alt_support -fih 5 -fiw 5 -rot 45 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeatureHistos.py -fn Figure_S3C -pt mutation -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f bestWTPeptideCount_I -fih 5 -fiw 5 -rot 45 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeatureHistos.py -fn Figure_S3E -pt mutation -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f bestWTMatchScore_I -fih 5 -fiw 5 -rot 45 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeatureHistos.py -fn Figure_S3G -pt mutation -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f CSCAPE_score -fih 5 -fiw 5 -rot 0 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeatureHistos.py -fn Figure_S3H -pt mutation -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f mutation_driver_statement_Intogen -fih 5 -fiw 5 -rot 60 -nt -ena"
echo $cmd
eval $cmd
cmd="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotFeatureHistos.py -fn Figure_S3I -pt mutation -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f nb_same_mutation_Intogen -fih 5 -fiw 5 -rot 45 -nt -log"
echo $cmd
eval $cmd

