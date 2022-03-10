import os.path

print("Appending previous log file Classifier_SVM_log.txt to Classifier_SVM_log_all.txt")

if os.path.isfile('/home/localadmin/Priorization/Logs/Classifier_SVM_log.txt'):
    with open('/home/localadmin/Priorization/Logs/Classifier_SVM_log.txt', "r") as last_out_file:
        with open('/home/localadmin/Priorization/Logs/Classifier_SVM_log_all.txt', "a+") as all_out_file:
            all_out_file.write(last_out_file.read())

