from sklearn.ensemble import VotingClassifier
from sklearn.metrics import make_scorer
from Classify_functions import *
import warnings
from sklearn.exceptions import UndefinedMetricWarning
from sklearn.preprocessing import *
import random

warnings.filterwarnings(action='ignore', category=UndefinedMetricWarning)

create_Rosenberg_plots = True
create_Tesla_plots = True
create_0YM1_plots = True
create_058C_plots = True
create_optim_plots = False
#normalizer = StandardScaler()
normalizer = QuantileTransformer()
#normalizer = None
verbose = 0
nr_clfs = 1
nr_iter = 10
optimize_function = optimize_svm_rdn
use_decision_score = True
#optimize_function = optimize_rf_rdn
#use_decision_score = False
#optimize_function = optimize_cart_rdn
#use_decision_score = False


random.seed(13071963)

#top100_scorer = make_scorer(nr_correct_top100, needs_threshold = True)
top100_scorer = make_scorer(sum_prob_correct, needs_threshold = True)

rosenberg_res_file = \
    '~/Documents/data/Prediction/Plots_Tesla_Rosenberg/Rosenberg_SVM_results.txt'

tesla1_res_file = \
    '~/Documents/data/Prediction/Plots_Tesla_Rosenberg/Tesla1_SVM_results.txt'

tesla3_res_file = \
    '~/Documents/data/Prediction/Plots_Tesla_Rosenberg/Tesla3_SVM_results.txt'

OYM1_res_file = \
    '~/Documents/data/Prediction/Plots_Tesla_Rosenberg/0YM1_SVM_results.txt'

O58C_res_file = \
    '~/Documents/data/Prediction/Plots_Tesla_Rosenberg/058C_SVM_results.txt'

df_R, y_R, df_T1, y_T1, df_T3, y_T3, df_0YM1, y_0YM1, df_058C, y_058C = load_all()

patients = np.array(df_R.patient)

cl_vars = ['rnaseq_ref_support','Sample_Tissue_expression_GTEx','MIN_MUT_RANK_CI_MIXMHC','WT_BEST_RANK_CI_MIXMHC',
           'COUNT_MUT_RANK_CI_MIXMHC','WT_RANK_CI_MIXMHC','MIN_MUT_RANK_CI_PRIME','WT_BEST_RANK_CI_PRIME',
           'COUNT_MUT_RANK_CI_PRIME','WT_RANK_CI_PRIME','CSCAPE_score','GTEx_all_tissues_expression_median',
           'bestWTMatchScore_I','bestWTMatchOverlap_I','bestWTPeptideCount_I','bestWT_Cancer_ProteinCnt_I',
           'rnaseq_gene_expression_quartile']

cl_vars = ['rnaseq_ref_support','Sample_Tissue_expression_GTEx','MIN_MUT_RANK_CI_MIXMHC','WT_BEST_RANK_CI_MIXMHC',
           'MIN_MUT_RANK_CI_PRIME','WT_BEST_RANK_CI_PRIME',
           'CSCAPE_score','GTEx_all_tissues_expression_median',
           'bestWTMatchScore_I','bestWTPeptideCount_I','bestWT_Cancer_ProteinCnt_I',
           'rnaseq_gene_expression_quartile']

# cl_vars = ['rnaseq_ref_support','Sample_Tissue_expression_GTEx','MIN_MUT_RANK_CI_MIXMHC','WT_BEST_RANK_CI_MIXMHC',
#            'COUNT_MUT_RANK_CI_MIXMHC','WT_RANK_CI_MIXMHC','MIN_MUT_RANK_CI_PRIME','WT_BEST_RANK_CI_PRIME',
#            'COUNT_MUT_RANK_CI_PRIME','WT_RANK_CI_PRIME','CSCAPE_score','GTEx_all_tissues_expression_median',
#            'rnaseq_gene_expression_quartile']

# cl_vars = ['rnaseq_ref_support','Sample_Tissue_expression_GTEx','MIN_MUT_RANK_CI_MIXMHC','WT_BEST_RANK_CI_MIXMHC',
#            'CSCAPE_score','GTEx_all_tissues_expression_median','bestWTMatchScore_I',
#            'bestWTPeptideCount_I','rnaseq_gene_expression_quartile']

# cl_vars = ['rnaseq_ref_support','Sample_Tissue_expression_GTEx','MIN_MUT_RANK_CI_MIXMHC','MIN_MUT_RANK_CI_PRIME',
#            'CSCAPE_score','GTEx_all_tissues_expression_median','bestWTMatchScore_I','bestWTMatchOverlap_I',
#            'bestWTPeptideCount_I','bestWT_Cancer_ProteinCnt_I','rnaseq_gene_expression_quartile']

print('Rosenberg data class counter: ' + str(Counter(y_R)))
if normalizer is not None:
    X_R = normalizer.fit_transform(df_R[cl_vars])
    X_T1 = normalizer.fit_transform(df_T1[cl_vars])
    X_T3 = normalizer.fit_transform(df_T3[cl_vars])
    X_0YM1 = normalizer.fit_transform(df_0YM1[cl_vars])
    X_058C = normalizer.fit_transform(df_058C[cl_vars])
else:
    X_R = df_R.loc[:,cl_vars].to_numpy()
    X_T1 = df_T1.loc[:,cl_vars].to_numpy()
    X_T3 = df_T3.loc[:,cl_vars].to_numpy()
    X_0YM1 = df_0YM1.loc[:,cl_vars].to_numpy()
    X_058C = df_058C.loc[:,cl_vars].to_numpy()

nr_correct_tot = 0
nr_immuno_tot = 0
nr_missed_tot = 0

y_pred_all = np.array([])
score_pred_all = np.array([])
cv_layer = np.array([])
y_all = np.array([])
X_all = None
p_all = np.array([])
id_all = np.array([])

splits = train_test_split(patients, y_R, 0.2)
cv_cnt = 0
cv_clfs = []
cv_weights = []
for train_idx, test_idx in splits:

    cv_cnt += 1

    X_train = X_R[train_idx,:]
    X_test = X_R[test_idx,:]

    y_train = np.array(df_R.loc[train_idx,'response'])
    y_test = np.array(df_R.loc[test_idx,'response'])
    p_test = np.array(df_R.loc[test_idx,'patient'])
    id_test = np.array(df_R.loc[test_idx,'peptide_id'])

    if np.sum(y_test == 1) == 0:
        continue

    best_clfs, best_scores, best_std = \
        optimize_function(X_train, y_train, top100_scorer, create_optim_plots, nr_clfs=nr_clfs, verbose=0,
                         nr_iter=nr_iter)

    estimators = []
    weights = []
    dec_fct = np.full(len(y_test),0.0)
    for i in range(len(best_clfs)):
        entry = ('svm_%d_%d' %(cv_cnt, i), best_clfs[i])
        estimators.append(entry)
        cv_clfs.append(entry)
        w = best_scores[i]-best_std[i]
        weights.append(w)
        cv_weights.append(w)
        best_clfs[i].fit(X_train, y_train)
        if use_decision_score: dec_fct += best_clfs[i].decision_function(X_test)

    voting_clf = VotingClassifier(estimators = estimators, voting='soft',weights=None)
    voting_clf.fit(X_train, y_train)

    y_all = np.append(y_all,y_test)
    if X_all is None:
        X_all = X_test
    else:
        X_all = np.concatenate((X_all, X_test))

    y_pred_all = np.append(y_pred_all, voting_clf.predict_proba(X_test)[:,1])
    score_pred_all = np.append(score_pred_all, dec_fct)

    cv_layer = np.append(cv_layer, np.full(len(p_test),cv_cnt))
    p_all = np.append(p_all, p_test)
    id_all = np.append(id_all, id_test)

    y_pred, nr_correct, nr_immuno, nr_missed = \
        test_classifier_on_patients(voting_clf, p_test, X_test, y_test, 20, verbose=1, prob= not use_decision_score)

    nr_correct_tot += nr_correct
    nr_immuno_tot += nr_immuno
    nr_missed_tot +=  nr_missed

    print('True positive cnt: %d, immunogenic cnt: %d, missed patient cnt: %d' %(nr_correct, nr_immuno, nr_missed))


print('True positive cnt: %d, immunogenic cnt: %d, missed patient cnt: %d' % (nr_correct_tot, nr_immuno_tot,
                                                                              nr_missed_tot))

if create_Rosenberg_plots:

    df_test = pd.DataFrame(data=X_all, columns=cl_vars)
    df_test.insert(len(cl_vars), "score", y_pred_all)
    df_test.insert(len(cl_vars) + 1, "response", y_all)
    df_test.insert(len(cl_vars) + 2, "patient", p_all)
    df_test.insert(len(cl_vars) + 3, "peptide_id", id_all)
    df_test.insert(len(cl_vars) + 4, "probability", score_pred_all)

    df_test.to_csv(rosenberg_res_file, sep = "\t", header=True, index=False)

    plot_class_results(df_test,"Rosenberg")
    plt.show(block=False)


voting_clf = VotingClassifier(estimators=cv_clfs, voting='soft', weights=cv_weights)
voting_clf.fit(X_R, y_R)

y_pred_T1, nr_correct, nr_immuno, nr_missed = \
    test_classifier_on_patients(voting_clf, np.full(df_T1.shape[0],'TESLA1'), X_T1, y_T1, 20, verbose=1, prob=True)

y_pred_T3, nr_correct, nr_immuno, nr_missed = \
    test_classifier_on_patients(voting_clf, np.full(df_T3.shape[0],'TESLA3'), X_T3, y_T3, 20, verbose=1, prob=True)

y_pred_0YM1, nr_correct, nr_immuno, nr_missed = \
    test_classifier_on_patients(voting_clf, np.full(df_0YM1.shape[0],'0YM1'), X_0YM1, y_0YM1, 20, verbose=1, prob=True)

y_pred_058C, nr_correct, nr_immuno, nr_missed = \
    test_classifier_on_patients(voting_clf, np.full(df_058C.shape[0], '058C'), X_058C, y_058C, 20, verbose=1, prob=True)

# if create_Rosenberg_plots:
#
#     df = df_R[cl_vars]
#     df.insert(len(cl_vars), "score", y_pred_all)
#     df.insert(len(cl_vars) + 1, "response", y_all)
#     df.insert(len(cl_vars) + 2, "patient", p_all)
#     df.insert(len(cl_vars) + 3, "peptide_id", id_all)
#     df.insert(len(cl_vars) + 4, "probability", score_pred_all)
#
#     df.to_csv(rosenberg_res_file, sep = "\t", header=True, index=False)
#
#     plot_class_results(df)
#     plt.show(block=False)


if create_Tesla_plots:

    df = pd.DataFrame(data=df_T1[cl_vars], columns=cl_vars)
    df.insert(len(cl_vars), "score", y_pred_T1)
    df.insert(len(cl_vars) + 1, "response", y_T1)
    df.insert(len(cl_vars) + 2, "patient", np.full(df.shape[0],'TESLA1'))
    df.insert(len(cl_vars) + 3, "peptide_id", np.full(df.shape[0],'TESLA1'))
    df.insert(len(cl_vars) + 4, "probability", y_pred_T1)

    df.to_csv(tesla1_res_file, sep = "\t", header=True, index=False)

    plot_class_results(df,"TESLA1")
    plt.show(block=False)

    df = pd.DataFrame(data=df_T3[cl_vars], columns=cl_vars)
    df.insert(len(cl_vars), "score", y_pred_T3)
    df.insert(len(cl_vars) + 1, "response", y_T3)
    df.insert(len(cl_vars) + 2, "patient", np.full(df.shape[0],'TESLA3'))
    df.insert(len(cl_vars) + 3, "peptide_id", np.full(df.shape[0],'TESLA3'))
    df.insert(len(cl_vars) + 4, "probability", y_pred_T3)

    df.to_csv(tesla3_res_file, sep = "\t", header=True, index=False)

    plot_class_results(df,"TESLA3")

    plt.show(block=False)


if create_0YM1_plots:

    df = pd.DataFrame(data=df_0YM1[cl_vars], columns=cl_vars)
    df.insert(len(cl_vars), "score", y_pred_0YM1)
    df.insert(len(cl_vars) + 1, "response", y_0YM1)
    df.insert(len(cl_vars) + 2, "patient", np.full(df.shape[0],'0YM1'))
    df.insert(len(cl_vars) + 3, "peptide_id", np.full(df.shape[0],'0YM1'))
    df.insert(len(cl_vars) + 4, "probability", y_pred_0YM1)

    df.to_csv(OYM1_res_file, sep = "\t", header=True, index=False)

    plot_class_results(df,"0MY1")
    plt.show(block=False)


if create_058C_plots:

    df = pd.DataFrame(data=df_058C[cl_vars], columns=cl_vars)
    df.insert(len(cl_vars), "score", y_pred_058C)
    df.insert(len(cl_vars) + 1, "response", y_058C)
    df.insert(len(cl_vars) + 2, "patient", np.full(df.shape[0], '058C'))
    df.insert(len(cl_vars) + 3, "peptide_id", np.full(df.shape[0], '058C'))
    df.insert(len(cl_vars) + 4, "probability", y_pred_058C)

    df.to_csv(O58C_res_file, sep="\t", header=True, index=False)

    plot_class_results(df,"058C")
    plt.show(block=False)


if create_Rosenberg_plots or create_Tesla_plots or create_058C_plots or create_0YM1_plots:
    plt.show()


# gamma : if overfitting -> decrease gamma. kernel(x','s) = exp[-gamma*(x-s)^2] : the higher gamma the more local the kernel
# C :if overfitting -> decrease C. min 1/2*w*w^t + C*sum_i(eta_i) : eta_i slack variable = (0 if class true',' dist to margin otherwise)
# class_weight : eg. class_weight={0: 1',' 1: 2}. C is replaced by C*class_weight[i]

# best so far {'C': 1.5',' 'class_weight': {1: 9}',' 'gamma': 1}
# Best params: {'C': 0.2750980053342904',' 'class_weight': {1: 11}',' 'gamma': 2.106785011385351',' 'kernel': 'rbf'}',' best precision: 0.11773887232449073
# Best params: {'C': 0.23091500814262222',' 'class_weight': {1: 10}',' 'gamma': 0.66825290705259',' 'kernel': 'rbf'}',' best precision: 0.2000549224220788
# 2885  C:0.213775       Gamma:0.73727      cl_weights:{1: 10}                      mean_precision:0.200055           0.400041                   1
# Best params: {'C': 0.06175348288740734',' 'class_weight': {1: 10}',' 'gamma': 3.609993861334124',' 'kernel': 'rbf'}',' nr_correct_top20: 2.0



