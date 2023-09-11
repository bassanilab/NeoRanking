source configure.sh

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotClassifierResults.py -d sampling_nr -fn Figure_S4A -re LR_train_all_Gartner_train_200_10000_short_*test.txt LR_train_all_Gartner_train_200_25000_short_*test.txt LR_train_all_Gartner_train_200_50000_short_*test.txt LR_train_all_Gartner_train_200_75000_short_*test.txt LR_train_all_Gartner_train_200_100000_short_09*test.txt LR_train_all_Gartner_train_200_150000_short_09*test.txt LR_train_all_Gartner_train_200_200000_short_09*test.txt -pt neopep -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'10\'000',1:'25\'000',2:'50\'000',3:'75\'000',4:'100\'000',5:'150\'000',6:'200\'000'}\" -a 0.02 -las 30 -tis 20 -rot 30 -tts 25 -ttp \"LR, \" -ds NCI_train -plt rank_score -oc cornflowerblue"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotClassifierResults.py -d hyperopt_iter -fn Figure_S4B -re LR_train_all_Gartner_train_10_*test.txt LR_train_all_Gartner_train_25_*test.txt LR_train_all_Gartner_train_50_*test.txt LR_train_all_Gartner_train_100_*test.txt LR_train_all_Gartner_train_150_*test.txt LR_train_all_Gartner_train_200_*test.txt LR_train_all_Gartner_train_250_*test.txt -pt neopep -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'10 iter',1:'25 iter',2:'50 iter',3:'100 iter',4:'150 iter',5:'200 iter',6:'250 iter'}\" -a 0.02 -las 30 -tis 20 -rot 30 -tts 25 -ttp \"LR, \" -ds NCI_train -plt rank_score -oc cornflowerblue"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotClassifierResults.py -d normalization -fn Figure_S4C -re LR_train_all_Gartner_train_200_100000_0.005_short*test.txt LR_train_all_Gartner_train_200_100000_0.005_p_short_*test.txt LR_train_all_Gartner_train_200_100000_0.005_z_short_*test.txt LR_train_all_Gartner_train_200_100000_0.005_i_short_*test.txt LR_train_all_Gartner_train_200_100000_0.005_n_short_*test.txt -pt neopep -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'quantile',1:'power',2:'z-score',3:'min-max',4:'none'}\" -a 0.02 -las 30 -tis 20 -rot 30 -tts 25 -ttp \"LR, \" -ds NCI_train -plt rank_score -oc cornflowerblue"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotClassifierResults.py -d normalization -fn Figure_S4D -re LR_train_all_Gartner_train_200_100000_0.005_short*test.txt LR_train_all_Gartner_train_200_100000_0.005_p_short_*test.txt LR_train_all_Gartner_train_200_100000_0.005_z_short_*test.txt LR_train_all_Gartner_train_200_100000_0.005_i_short_*test.txt LR_train_all_Gartner_train_200_100000_0.005_n_short_*test.txt -pt neopep -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'quantile',1:'power',2:'z-score',3:'min-max',4:'none'}\" -a 0.02 -las 30 -tis 20 -rot 30 -tts 25 -ttp \"LR, \" -ds NCI_train -plt topn_counts -oc cornflowerblue"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotClassifierResults.py -d classifier_comp -fn Figure_S4E -re LR_*test.txt SVM_*test.txt SVM-lin_*test.txt XGBoost_*test.txt CatBoost_*test.txt -pt neopep -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR',1:'SVM-RBF',2:'SVM-Linear',3:'XGBoost',4:'CatBoost'}\" -a 0.02 -las 30 -tis 20 -rot 30 -tts 25 -ds NCI_train -plt rank_score -oc cornflowerblue"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotClassifierResults.py -d classifier_comp -fn Figure_S4F -re LR_*test.txt SVM_*test.txt SVM-lin_*test.txt XGBoost_*test.txt CatBoost_*test.txt -pt neopep -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR',1:'SVM-RBF',2:'SVM-Linear',3:'XGBoost',4:'CatBoost'}\" -a 0.02 -las 30 -tis 20 -rot 30 -tts 25 -ds NCI_train -plt topn_counts"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotClassifierPCA.py -d classifier_comp -fn Figure_S4G -re LR_*test.txt SVM_*test.txt SVM-lin_*test.txt XGBoost_*test.txt CatBoost_*test.txt -pt neopep -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR',1:'SVM-RBF',2:'SVM-Linear',3:'XGBoost',4:'CatBoost'}\" -a 0.02 -las 30 -tis 20 -rot 30 -tts 25 -ds NCI_train -cm \"{'LR':0,'SVM-RBF':1,'SVM-Linear':2,'XGBoost':3,'CatBoost':4}\""
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotClassifierResults.py -d hitide_train -fn Figure_S4H -re LR_train_all_Gartner*_test.txt LR_train_all_HiTIDE*_test.txt -pt neopep -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR NCI-train',1:'LR HiTIDE'}\" -a 0.02 -las 30 -tis 20 -rot 30 -tts 25 -ds NCI_test -plt topn_counts"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotClassifierResults.py -d hitide_train -fn Figure_S4I -re LR_train_all_Gartner*_test.txt LR_train_all_HiTIDE*_test.txt -pt neopep -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR NCI-train',1:'LR HiTIDE'}\" -a 0.02 -las 30 -tis 20 -rot 30 -tts 25 -ds TESLA -plt topn_counts"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotClassifierResults.py -d hitide_train -fn Figure_S4J -re LR_train_all_Gartner*_test.txt LR_train_all_HiTIDE*_test.txt -pt neopep -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR NCI-train',1:'LR HiTIDE'}\" -a 0.02 -las 30 -tis 20 -rot 30 -tts 25 -ds HiTIDE -plt topn_counts"
echo $CMD
eval $CMD

#CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotShapleyImportance.py -d hitide_train -c1 LR_*Gartner*.sav -c2 LR_*HiTIDE*.sav -ds1 NCI_train -ds2 HiTIDE -pt neopep -ft pdf -fn Figure_S4K -las 30 -tis 23 -tts 25 -les 20"
#echo $CMD
#eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotClassifierResults.py -d ipmsdb_intogen_contribution -fn Figure_S4L -re LR_train_all*_short_*test.txt LR_train_noMS_*_short_*test.txt LR_train_noINT*_short_*test.txt LR_train_noMSINT*_short_*test.txt -pt neopep -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'ALL',1:'No ipMSDB',2:'No Intogen',3:'No ipMSDB & Intogen'}\" -a 0.02 -rot 75 -las 30 -tis 18 -les 15 -sr -ds NCI_test -plt rank_score -oc cornflowerblue"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotClassifierResults.py -fn Figure_S4N -re LR_paper_data_*_mutation_*test.txt XGBoost_paper_data_*_mutation_*test.txt Voting_classifier_mutation_0.50_test.txt -pt mutation -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR',1:'XGBoost',2:'Voting'}\" -a 0.02 -rot 30 -las 30 -tis 18 -les 15 -sr -o \"LR,XGBoost,Voting\" -ds TESLA -plt topn_counts"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotClassifierResults.py -fn Figure_S4M -re LR_paper_data_*_mutation_*test.txt XGBoost_paper_data_*_mutation_*test.txt Voting_classifier_mutation_0.50_test.txt -pt mutation -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR',1:'XGBoost',2:'Voting'}\" -a 0.02 -rot 30 -las 30 -tis 18 -les 15 -sr -o \"LR,XGBoost,Voting\" -ds HiTIDE -plt topn_counts"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotClassifierResults.py -d ipmsdb_intogen_contribution -fn Figure_S4O -re LR_train_all*_long_*test.txt LR_train_noMS_*_long_*test.txt LR_train_noINT*_long_*test.txt LR_train_noMSINT*_long_*test.txt -pt neopep -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'ALL',1:'No ipMSDB',2:'No Intogen',3:'No ipMSDB & Intogen'}\" -a 0.02 -rot 75 -las 30 -tis 18 -les 15 -sr -ds NCI_test -plt rank_score -oc cornflowerblue"
echo $CMD
eval $CMD

