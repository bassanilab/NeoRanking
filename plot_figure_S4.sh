cmd="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 Scripts/paper/PlotFeatureHistos.py -fn Suppl_Figure_4A -pt mutation -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f rnaseq_TPM -fih 5 -fiw 5 -rot 45 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 Scripts/paper/PlotFeatureHistos.py -fn Suppl_Figure_4B -pt mutation -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f rnaseq_alt_support -fih 5 -fiw 5 -rot 45 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 Scripts/paper/PlotFeatureHistos.py -fn Suppl_Figure_4C -pt mutation -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f bestWTPeptideCount_I -fih 5 -fiw 5 -rot 45 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 Scripts/paper/PlotFeatureHistos.py -fn Suppl_Figure_4E -pt mutation -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f bestWTMatchScore_I -fih 5 -fiw 5 -rot 45 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 Scripts/paper/PlotFeatureHistos.py -fn Suppl_Figure_4G -pt mutation -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f CSCAPE_score -fih 5 -fiw 5 -rot 0 -nt"
echo $cmd
eval $cmd
cmd="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 Scripts/paper/PlotFeatureHistos.py -fn Suppl_Figure_4H -pt mutation -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f mutation_driver_statement_Intogen -fih 5 -fiw 5 -rot 60 -nt -ena"
echo $cmd
eval $cmd
cmd="PYTHONPATH=/home/localadmin/Priorization/paper_code python3 Scripts/paper/PlotFeatureHistos.py -fn Suppl_Figure_4I -pt mutation -ds NCI -ds TESLA -ds HiTIDE -ft pdf -f nb_same_mutation_Intogen -fih 5 -fiw 5 -rot 45 -nt -log"
echo $cmd
eval $cmd

