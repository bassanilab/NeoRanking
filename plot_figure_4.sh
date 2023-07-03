source configure.sh

# plot Gartner comparison
CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotMultiClassifierOutput.py -fn Figure_4A -re LR_train_all_*_long_*test.txt XGBoost_*_long_*test.txt Voting_classifier_0.50_long_test.txt -pt mutation -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR',1:'XGBoost',2:'Voting'}\" -a 0.02 -rot 30 -las 30 -tis 18 -les 15 -sr -ga -o \"Gartner et al.,LR,XGBoost,Voting\" -ds NCI_test -plt topn_counts"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotMultiClassifierOutput.py -fn Figure_4B -re LR_train_all_*_long_*test.txt XGBoost_*_long_*test.txt Voting_classifier_0.50_long_test.txt -pt mutation -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR',1:'XGBoost',2:'Voting'}\" -a 0.02 -rot 30 -las 30 -tis 18 -les 15 -sr -o \"LR,XGBoost,Voting\" -ds TESLA -plt topn_counts"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotMultiClassifierOutput.py -fn Figure_4C -re LR_train_all_*_long_*test.txt XGBoost_*_long_*test.txt Voting_classifier_0.50_long_test.txt -pt mutation -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'LR',1:'XGBoost',2:'Voting'}\" -a 0.02 -rot 30 -las 30 -tis 18 -les 15 -sr -o \"LR,XGBoost,Voting\" -ds HiTIDE -plt topn_counts"
echo $CMD
eval $CMD

CMD="PYTHONPATH=$NEORANKING_CODE python3 Scripts/paper/PlotMultiClassifierOutput.py -fn Figure_4E -re LR_train_all*_long_*test.txt LR_train_noMS*_long_*test.txt LR_train_noINT*_long_*test.txt LR_train_noMSINT*_long_*test.txt -pt mutation -g \"{'TESLA':'TESLA','HiTIDE':'HiTIDE','NCI_test':'NCI-test','NCI_train':'NCI-train'}\" -cln \"{0:'ALL',1:'No ipMSDB',2:'No Intogen',3:'No ipMSDB & Intogen'}\" -a 0.02 -rot 75 -las 30 -tis 18 -les 15 -sr -ds NCI_test -plt rank_score -oc cornflowerblue"
echo $CMD
eval $CMD

