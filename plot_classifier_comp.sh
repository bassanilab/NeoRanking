## short

### plot classifier comparison
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers -fp Clf_transfer_HiTIDE_NCI -re LR_*_Gartner_*test.txt Voting_classifier_1.00_test.txt Voting_classifier_0.50_test.txt -pt short -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR NCI-train',1:'LR HiTIDE',2:'Voting w=0.5'}\" -a 0.02 -las 30 -tis 18 -les 15"
#echo $CMD
#eval $CMD

### plot classifier comparison
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_eval_transferL -fp Clf_HiTIDE_NCI_comp -re LR_train_all_Gartner*_test.txt LR_train_all_HiTIDE*_test.txt -pt short -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR NCI-train',1:'LR HiTIDE'}\" -a 0.02 -las 30 -tis 18 -les 15 -o \"LR NCI-train,LR HiTIDE\""
echo $CMD
eval $CMD

### plot classifier comparison
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_eval_transferL -fp Clf_transfer_HiTIDE_NCI_NeoDisc -re Voting_classifier_0.00_test.txt Voting_classifier_0.40_test.txt Voting_classifier_1.00_test.txt -pt short -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'voting weight=0.0',1:'voting weight=0.4',2:'voting weight=1.0'}\" -a 0.02 -rot 60 -las 30 -tis 18 -les 15 -nd"
#echo $CMD
#eval $CMD

### plot classifier comparison
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_eval_transferL -fp Clf_transfer_HiTIDE_NCI_weights -re Voting_classifier_0.00_test.txt Voting_classifier_0.10_test.txt Voting_classifier_0.20_test.txt Voting_classifier_0.30_test.txt Voting_classifier_0.40_test.txt Voting_classifier_0.50_test.txt Voting_classifier_0.60_test.txt Voting_classifier_0.70_test.txt Voting_classifier_0.80_test.txt Voting_classifier_0.90_test.txt Voting_classifier_1.00_test.txt -pt short -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'0.0',1:'0.1',2:'0.2',3:'0.3',4:'0.4',5:'0.5',6:'0.6',7:'0.7',8:'0.8',9:'0.9',10:'1.0'}\" -a 0.02 -las 30 -tis 18 -les 15 -o \"0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0\" -xl Weight -bar -oc cornflowerblue -ylim \"{'HiTIDE':10}\""
#echo $CMD
#eval $CMD

## plot classifier comparison for LR, SVM, XGBoost and CatBoost
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_eval_short -fp Clf_comp -re LR_*test.txt  SVM_*test.txt SVM-lin_*test.txt XGBoost_*test.txt CatBoost_*test.txt -pt short -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR',1:'SVM-RBF',2:'SVM-Linear',3:'XGBoost',4:'CatBoost'}\" -a 0.02 -las 30 -tis 18 -les 15 -hy -ylim \"{'NCI_train':25}\" "
#echo $CMD
#eval $CMD

## plot classifier comparison with mixmhc and netmhc ranking, LR, SVM, XGBoost, CatBoost and voting
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_eval_short -fp Clf_comp_mixmhc_voting -re LR_*test.txt  SVM_*test.txt SVM-lin_*test.txt XGBoost_*test.txt CatBoost_*test.txt Voting_classifier_test.txt SimpleRanking_short_mixmhc_test.txt SimpleRanking_short_netmhc_test.txt -pt short -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR',1:'SVM-RBF',2:'SVM-Linear',3:'XGBoost',4:'CatBoost',5:'Voting',6:'MixMHC+RNA',7:'NetMHC+RNA'}\" -a 0.02 -rot 75 -las 30 -tis 18 -les 15 -hy -sr -o \"MixMHC+RNA,NetMHC+RNA,LR,SVM-RBF,SVM-Linear,XGBoost,CatBoost,Voting\""
#echo $CMD
#eval $CMD

## plot classifier comparison with mixmhc and netmhc ranking, LR, SVM, XGBoost, CatBoost and voting
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_eval_short -fp Clf_comp_mixmhc_voting -re LR_*test.txt  XGBoost_*test.txt Voting_classifier_0.50_test.txt SimpleRanking_short_mixmhc_test.txt SimpleRanking_short_netmhc_test.txt -pt short -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR',1:'XGBoost',2:'Voting',3:'MixMHC+RNA',4:'NetMHC+RNA'}\" -a 0.02 -rot 75 -las 30 -tis 18 -les 15 -hy -sr -ga -o \"MixMHC+RNA,NetMHC+RNA,Gartner et al.\nranking,LR,XGBoost,Voting\""
#echo $CMD
#eval $CMD

# plot classifier comparison with simple ranking and voting
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_eval_short -fp Clf_comp_sira_voting -re LR_*test.txt XGBoost_*test.txt Voting*_test.txt SimpleRanking_short_affinity_rank_test.txt -pt short -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR',1:'XGBoost',2:'Voting',3:'MixMHC+RNA'}\" -a 0.02 -rot 40 -las 30 -tis 18 -les 15 -sr -o \"MixMHC+RNA,LR,XGBoost,Voting\""
#echo $CMD
#eval $CMD


# plot classifier comparison with simple ranking and voting
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_eval-short -fp Clf_comp_filter_mut_LR_voting -re Voting_classifier_test.txt Voting_classifier_rnaseq_test.txt Voting_classifier_filter_0.50_50_5_test.txt -pt short -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'Voting',1:'Voting+FilterRNA',2:'Voting+FilterMut'}\" -a 0.02 -rot 40 -las 30 -tis 18 -les 15 -sr -o \"Voting,Voting+FilterRNA,Voting+FilterMut\""
#echo $CMD
#eval $CMD

## plot classifier comparison with simple ranking and voting
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_eval_short -fp Clf_comp_filter_ipMSDB_voting -re LR_*test.txt XGBoost_*test.txt Voting_classifier_test.txt Voting_classifier_0.50_test.txt -pt short -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR',1:'XGBoost',2:'Voting',3:'Voting+ipMSDB'}\" -a 0.02 -rot 40 -las 30 -tis 18 -les 15 -sr -o \"LR,XGBoost,Voting,Voting+ipMSDB\""
#echo $CMD
#eval $CMD

# plot feature comparison
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_feature -fp LR_fsets_comp -re LR_train_set1_*test.txt LR_train_set2_*test.txt LR_train_set3_*test.txt LR_train_set4_*test.txt LR_train_set5_*test.txt LR_train_set6*_test.txt LR_train_set7*_test.txt LR_train_set8*_test.txt LR_train_all*_test.txt LR_train_set9*_test.txt LR_train_set10*_test.txt LR_train_set11*_test.txt LR_train_set12*_test.txt LR_train_set13*_test.txt LR_train_set1_Gartner_train_50_5000_*_test.txt -pt short -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'MixMHC',1:'NetMHC',2:'RNA exp',3:'TCGA exp',4:'RNA cov',5:'NetStab/TAP/Chop',6:'ipMSDB',7:'Intogen',8:'all',9:'RNA+ipMSDB',10:'RNA+Intogen',11:'RNA+Intogen',12:'Intogen+MatchType',13:'MatchType',14:'MixMHC 5000'}\" -a 0.02 -las 20 -tis 18 -les 15 -hy -rot 40"
#echo $CMD
#eval $CMD

# plot LR performance with variable nr negatives for LR
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_eval_nn -fp LR_sampling -re LR_train_all_Gartner_train_200_10000_short_*test.txt LR_train_all_Gartner_train_200_25000_short_*test.txt LR_train_all_Gartner_train_200_50000_short_*test.txt LR_train_all_Gartner_train_200_75000_short_*test.txt LR_train_all_Gartner_train_200_100000_short_09*test.txt LR_train_all_Gartner_train_200_150000_short_09*test.txt LR_train_all_Gartner_train_200_200000_short_09*test.txt -pt short -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'10\'000',1:'25\'000',2:'50\'000',3:'75\'000',4:'100\'000',5:'150\'000',6:'200\'000'}\" -a 0.02 -las 30 -tis 20 -rot 30 -tts 25 -ttp \"LR, \" -hy -oc cornflowerblue" 
#echo $CMD
#eval $CMD

# plot LR with variable alpha
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_alpha -fp LR_alpha -re LR_train_all_Gartner_train_200_100000_0.001_short_09.08.2022-12*test.txt LR_train_all_Gartner_train_200_100000_0.005_short_09.08.2022-17*test.txt LR_train_all_Gartner_train_200_100000_0.01_short_09.08.2022-14*test.txt -pt short -g  \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'alpha=0.001',1:'alpha=0.005',2:'alpha=0.01'}\" -a 0.02 -tis 15 -las 30 -rot 20 -tts 25 -ttp \"LR, \" -hy"
#echo $CMD
#eval $CMD

# plot LR with different normalization
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_norm -fp LR_norm -re LR_train_all_Gartner_train_200_100000_0.005_short*test.txt LR_train_all_Gartner_train_200_100000_0.005_p_short_*test.txt LR_train_all_Gartner_train_200_100000_0.005_z_short_*test.txt LR_train_all_Gartner_train_200_100000_0.005_i_short_*test.txt LR_train_all_Gartner_train_200_100000_0.005_n_short_*test.txt -pt short -g  \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'quantile',1:'power',2:'z-score',3:'min-max',4:'none'}\" -a 0.02 -tis 15 -hy -las 30 -tts 25 -ttp \"LR, \" -oc cornflowerblue -ylim \"{'NCI_train':25}\" "
#echo $CMD
#eval $CMD

# plot Gartner comparsion
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_eval_short -fp LR_XGB_Gartner_comp -re LR_*test.txt XGBoost_*test.txt Voting_classifier_0.50_test.txt Voting_classifier_filter_0.50_50_15_test.txt -pt short -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR',1:'XGBoost',2:'Voting',3:'Voting+FilterMut'}\" -a 0.02 -rot 45 -las 30 -tis 18 -les 15"
#echo $CMD
#eval $CMD 

# plot Gartner comparsion
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_eval_short -fp LR_XGB_Gartner_comp -re LR_*test.txt XGBoost_*test.txt Voting_classifier_0.50_test.txt Voting_classifier_filter_0.50_50_3_test.txt -pt short -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR',1:'XGBoost',2:'Voting',3:'Voting+FilterMut'}\" -a 0.02 -rot 45 -las 30 -tis 18 -les 15 -ga"
#echo $CMD
#eval $CMD 

# plot Gartner comparsion
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_eval_short -fp LR_XGB_Gartner_comp -re LR_*test.txt XGBoost_*test.txt Voting_classifier_0.50_test.txt -pt short -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR',1:'XGBoost',2:'Voting'}\" -a 0.02 -rot 45 -las 30 -tis 18 -les 15 -ga"
#echo $CMD
#eval $CMD 

# plot LR performance with variable nr hyperopt iterations
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_iter -fp LR_iter -re LR_train_all_Gartner_train_10_*test.txt LR_train_all_Gartner_train_25_*test.txt LR_train_all_Gartner_train_50_*test.txt LR_train_all_Gartner_train_100_*test.txt LR_train_all_Gartner_train_150_*test.txt LR_train_all_Gartner_train_200_*test.txt LR_train_all_Gartner_train_250_*test.txt -pt short -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'10 iter',1:'25 iter',2:'50 iter',3:'100 iter',4:'150 iter',5:'200 iter',6:'250 iter'}\" -a 0.02 -tis 15 -hy -las 30 -tts 25 -ttp \"LR, \" -oc cornflowerblue" 
#echo $CMD
#eval $CMD

## plot classifier comparison
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_fsel_short -fp Clf_feature_subset_short -re LR_train_all_**test.txt LR_train_noMS_*test.txt LR_train_noINT_*test.txt LR_train_noMSINT_*test.txt -pt short -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'ALL',1:'No ipMSDB',2:'No Intogen',3:'No ipMSDB & Intogen'}\" -a 0.02 -las 30 -tis 18 -les 15 -oc cornflowerblue"
#echo $CMD
#eval $CMD

## plot classifier comparison
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_eval_short /home/localadmin/Priorization/Classifiers_ipmsdb_qc_short /home/localadmin/Priorization/Classifiers_ipmsdb_noqc_short /home/localadmin/Priorization/Classifiers_eval_short /home/localadmin/Priorization/Classifiers_ipmsdb_qc_short -fp Clf_ipmsdb_qc_short -re LR_train_all_**test.txt LR_train_all_*test.txt LR_train_all_**test.txt XGBoost_train_all_*test.txt XGBoost_train_all_*test.txt -pt short -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR old ipMSDB',1:'LR new ipMSDB QC',2:'LR new ipMSDB no QC',3:'XGBoost old ipMSDB',4:'XGBoost new ipMSDB QC'}\" -a 0.02 -las 30 -tis 18 -les 15"
#echo $CMD
#eval $CMD

###############################################################################################################################################################################
## long

## plot classifier comparison
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_eval_long -fp Clf_comp_long -re LR_*test.txt XGBoost_*test.txt -pt long -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR',1:'XGBoost'}\" -a 0.02 -las 30 -tis 18 -les 15 -cm \"{'LR':0,'XGBoost':3}\" "
#echo $CMD
#eval $CMD

## plot classifier comparison
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_eval_long -fp Clf_comp_voting_long -re LR_*test.txt XGBoost_*test.txt Voting_classifier_0.50_test.txt -pt long -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR',1:'XGBoost',2:'Voting'}\" -a 0.02 -las 30 -tis 18 -les 15"
#echo $CMD
#eval $CMD

## plot classifier comparison
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_eval_long -fp Clf_comp_gartner_long -re LR_*test.txt XGBoost_*test.txt Voting_classifier_0.00_test.txt -pt long -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR',1:'XGBoost',2:'Voting'}\" -a 0.02 -las 30 -tis 18 -les 15 -ga"
#echo $CMD
#eval $CMD

## plot classifier comparison
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_fsel_long -fp Clf_feature_subset_long -re LR_train_all_**test.txt LR_train_noMS_*test.txt LR_train_noINT_*test.txt LR_train_noMSINT_*test.txt  -pt long -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'ALL',1:'No ipMSDB',2:'No Intogen',3:'No ipMSDB & Intogen'}\" -a 0.02 -las 30 -tis 18 -les 15 -oc cornflowerblue"
#echo $CMD
#eval $CMD

## plot classifier comparison
CMD="PYTHONPATH=/home/localadmin/Priorization/python python3 /home/localadmin/Priorization/python/Scripts/paper/PlotMultiClassifierOutput.py -d /home/localadmin/Priorization/Classifiers_eval_long /home/localadmin/Priorization/Classifiers_ipmsdb_qc_long /home/localadmin/Priorization/Classifiers_ipmsdb_noqc_long -fp Clf_ipmsdb_qc_long -re LR_train_all_**test.txt LR_train_all_*test.txt LR_train_all_**test.txt -pt long -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR old ipMSDB',1:'LR new ipMSDB QC',2:'LR new ipMSDB no QC'}\" -a 0.02 -las 30 -tis 18 -les 15"
#echo $CMD
#eval $CMD


