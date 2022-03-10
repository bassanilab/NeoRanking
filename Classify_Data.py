from sklearn.ensemble import VotingClassifier
from sklearn.metrics import make_scorer
from Classify_functions import *
import warnings
from sklearn.exceptions import UndefinedMetricWarning
from sklearn.preprocessing import *
import random
from UmapClusterProb import UmapClusterProb
import matplotlib.pyplot as plt
import pickle
from Utils.Parameters import *
from datetime import datetime

warnings.filterwarnings(action='ignore', category=UndefinedMetricWarning)

train_on_Rosenberg = False
create_Rosenberg_plots = False
create_Tesla_plots = False
create_0YM1_plots = True
create_058C_plots = False
create_optim_plots = False
#normalizer = StandardScaler()
normalizer = QuantileTransformer()
#normalizer = None
verbose = 1
nr_clfs = 1
nr_iter = 100
optimize_function = optimize_svm_rdn
#optimize_function = optimize_svmlin_rdn
#optimize_function = optimize_adaboost_rdn
use_decision_score = False
#optimize_function = optimize_rf_rdn
#optimize_function = optimize_cart_rdn
umaper = None


random.seed(13071963)

#top100_scorer = make_scorer(nr_correct_top100, needs_threshold = True)
top100_scorer = make_scorer(sum_rank_correct, needs_threshold = True)

class_out_file_tag = 'SVM_results'
parameters = Parameters()

df_R, y_R, df_T1, y_T1, df_T2, y_T2, df_T3, y_T3, df_T4, y_T4, df_T8, y_T8, df_T9, y_T9, \
df_T12, y_T12, df_T16, y_T16, df_0YM1, y_0YM1, df_058C, y_058C = \
    load_all('netmhc_stab_chop_mbp_tcr')

patients = np.array(df_R.patient)

cl_vars_1 = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx', 'MIN_MUT_RANK_CI_MIXMHC',
            'COUNT_MUT_RANK_CI_MIXMHC', 'WT_RANK_CI_MIXMHC', 'MIN_MUT_RANK_CI_PRIME',  'WT_BEST_RANK_CI_PRIME',
            'COUNT_MUT_RANK_CI_PRIME', 'WT_RANK_CI_PRIME', 'CSCAPE_score', 'GTEx_all_tissues_expression_median',
            'bestWTMatchScore_I', 'bestWTMatchOverlap_I', 'bestWTPeptideCount_I', 'bestWT_Cancer_ProteinCnt_I',
            'rnaseq_gene_expression_quartile', 'DAI']

cl_vars_2 = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx', 'RANK_CI_MIXMHC_DIFF',
           'MIN_MUT_RANK_CI_MIXMHC', 'COUNT_MUT_RANK_CI_MIXMHC', 'MIN_MUT_RANK_CI_PRIME',
           'CSCAPE_score', 'GTEx_all_tissues_expression_median',
           'bestWTMatchScore_I', 'bestWTPeptideCount_I', 'rnaseq_gene_expression_quartile',
           'mut_Rank_BA_0', 'mut_Rank_EL_0', 'mut_Rank_BA_1', 'mut_Rank_EL_1', 'mut_nr_weak_binding_alleles_0',
           'mut_nr_weak_binders_0']

cl_vars_3 = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx', 'RANK_CI_MIXMHC_DIFF',
           'MIN_MUT_RANK_CI_MIXMHC', 'COUNT_MUT_RANK_CI_MIXMHC', 'MIN_MUT_RANK_CI_PRIME',
           'CSCAPE_score', 'GTEx_all_tissues_expression_median',
           'bestWTMatchScore_I', 'bestWTPeptideCount_I', 'rnaseq_gene_expression_quartile',
           'mut_Rank_BA_0', 'mut_Rank_EL_0', 'next_best_BA_mut_ranks', 'next_best_EL_mut_ranks', 'mut_nr_weak_binding_alleles_0',
           'mut_nr_weak_binders_0']

cl_vars_4 = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx', 'RANK_CI_MIXMHC_DIFF',
           'MIN_MUT_RANK_CI_MIXMHC', 'COUNT_MUT_RANK_CI_MIXMHC', 'MIN_MUT_RANK_CI_PRIME',
           'CSCAPE_score', 'GTEx_all_tissues_expression_median',
           'bestWTMatchScore_I', 'bestWTPeptideCount_I', 'rnaseq_gene_expression_quartile',
           'mut_Rank_BA_0', 'mut_Rank_EL_0', 'mut_Rank_BA_1', 'mut_Rank_EL_1', 'mut_nr_weak_binding_alleles_0',
           'mut_nr_weak_binders_0', 'DAI', 'mut_Rank_Stab_0', 'mut_Rank_Stab_1']

cl_vars_5 = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx', 'RANK_CI_MIXMHC_DIFF',
           'MIN_MUT_RANK_CI_MIXMHC', 'COUNT_MUT_RANK_CI_MIXMHC', 'MIN_MUT_RANK_CI_PRIME',
           'CSCAPE_score', 'GTEx_all_tissues_expression_median',
           'bestWTMatchScore_I', 'bestWTPeptideCount_I', 'rnaseq_gene_expression_quartile',
           'mut_Rank_EL_0', 'mut_Rank_EL_1', 'next_best_EL_mut_ranks', 'mut_nr_weak_binding_alleles_0',
           'mut_nr_weak_binders_0']

cl_vars_6 = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx', 'RANK_CI_MIXMHC_DIFF',
           'MIN_MUT_RANK_CI_MIXMHC', 'COUNT_MUT_RANK_CI_MIXMHC', 'MIN_MUT_RANK_CI_PRIME',
           'CSCAPE_score', 'GTEx_all_tissues_expression_median',
           'bestWTMatchScore_I', 'bestWTPeptideCount_I', 'rnaseq_gene_expression_quartile',
           'mut_Rank_BA_0', 'mut_Rank_EL_0', 'mut_Rank_BA_1', 'mut_Rank_EL_1', 'mut_nr_weak_binding_alleles_0',
           'mut_nr_weak_binders_0', 'DAI', 'mut_Rank_Stab_0', 'mut_Rank_Stab_1', 'mut_netchop_Ct_score_0',
           'mut_netchop_Ct_score_1']


cl_vars_7 = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx', 'RANK_CI_MIXMHC_DIFF',
           'MIN_MUT_RANK_CI_MIXMHC', 'COUNT_MUT_RANK_CI_MIXMHC', 'MIN_MUT_RANK_CI_PRIME',
           'CSCAPE_score', 'GTEx_all_tissues_expression_median',
           'bestWTMatchScore_I', 'bestWTPeptideCount_I', 'rnaseq_gene_expression_quartile',
           'mut_Rank_BA_0', 'mut_Rank_EL_0', 'DAI']

cl_vars_8 = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx', 'RANK_CI_MIXMHC_DIFF',
           'MIN_MUT_RANK_CI_MIXMHC', 'COUNT_MUT_RANK_CI_MIXMHC', 'MIN_MUT_RANK_CI_PRIME',
           'CSCAPE_score', 'GTEx_all_tissues_expression_median',
           'bestWTMatchScore_I', 'bestWTPeptideCount_I', 'rnaseq_gene_expression_quartile',
           'mut_Rank_BA_0', 'mut_Rank_EL_0', 'mut_Rank_BA_1', 'mut_Rank_EL_1', 'mut_nr_weak_binding_alleles_0',
           'mut_nr_weak_binders_0']

cl_vars_9 = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx', 'MIN_MUT_RANK_CI_PRIME',
           'CSCAPE_score', 'GTEx_all_tissues_expression_median',
           'bestWTMatchScore_I', 'bestWTPeptideCount_I', 'rnaseq_gene_expression_quartile',
           'mut_Rank_BA_0', 'mut_Rank_EL_0', 'mut_Rank_BA_1', 'mut_Rank_EL_1', 'mut_nr_weak_binding_alleles_0',
           'mut_nr_weak_binders_0', 'DAI']

cl_vars_10 = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx', 'RANK_CI_MIXMHC_DIFF',
           'MIN_MUT_RANK_CI_MIXMHC', 'COUNT_MUT_RANK_CI_MIXMHC', 'MIN_MUT_RANK_CI_PRIME',
           'CSCAPE_score', 'GTEx_all_tissues_expression_median',
           'bestWTMatchScore_I', 'bestWTPeptideCount_I', 'rnaseq_gene_expression_quartile',
           'mut_Rank_EL_0', 'next_best_EL_mut_ranks', 'mut_nr_weak_binding_alleles_0',
           'mut_nr_weak_binders_0', 'DAI']

cl_vars_11 = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx', 'RANK_CI_MIXMHC_DIFF',
           'MIN_MUT_RANK_CI_MIXMHC', 'COUNT_MUT_RANK_CI_MIXMHC', 'MIN_MUT_RANK_CI_PRIME',
           'CSCAPE_score', 'GTEx_all_tissues_expression_median',
           'bestWTMatchScore_I', 'bestWTPeptideCount_I', 'rnaseq_gene_expression_quartile',
           'mut_Rank_BA_0', 'mut_Rank_EL_0', 'next_best_BA_mut_ranks', 'next_best_EL_mut_ranks', 'mut_nr_weak_binding_alleles_0',
           'mut_nr_weak_binders_0', 'DAI', 'mut_Rank_Stab_0', 'mut_Rank_Stab_1', 'mut_netchop_Ct_score_0',
           'mut_netchop_Ct_score_1']

cl_vars_12 = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx', 'MIN_MUT_RANK_CI_PRIME',
           'CSCAPE_score', 'GTEx_all_tissues_expression_median',
           'bestWTMatchScore_I', 'bestWTPeptideCount_I', 'rnaseq_gene_expression_quartile',
           'mut_Rank_BA_0', 'mut_Rank_EL_0', 'next_best_BA_mut_ranks', 'next_best_EL_mut_ranks', 'mut_nr_weak_binding_alleles_0',
           'mut_nr_weak_binders_0', 'DAI', 'mut_Rank_Stab_0', 'mut_Rank_Stab_1', 'mut_netchop_Ct_score_0',
           'mut_netchop_Ct_score_1']

cl_vars_13 = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx', 'MIN_MUT_RANK_CI_PRIME',
           'CSCAPE_score', 'GTEx_all_tissues_expression_median',
           'bestWTMatchScore_I', 'bestWTPeptideCount_I', 'rnaseq_gene_expression_quartile',
           'mut_Rank_BA_0', 'mut_Rank_EL_0', 'next_best_BA_mut_ranks', 'next_best_EL_mut_ranks', 'mut_nr_weak_binding_alleles_0',
           'mut_nr_weak_binders_0', 'mut_Rank_Stab_0', 'mut_Rank_Stab_1', 'mut_netchop_Ct_score_0']

cl_vars_14 = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx',
           'CSCAPE_score', 'GTEx_all_tissues_expression_median',
           'bestWTMatchScore_I', 'bestWTPeptideCount_I', 'rnaseq_gene_expression_quartile',
           'mut_Rank_BA_0', 'mut_Rank_EL_0', 'next_best_BA_mut_ranks', 'next_best_EL_mut_ranks', 'mut_nr_weak_binding_alleles_0',
           'mut_nr_weak_binders_0', 'mut_Rank_Stab_0', 'mut_Rank_Stab_1', 'mut_netchop_Ct_score_0']

cl_vars_15 = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx', 'MIN_MUT_RANK_CI_PRIME',
           'GTEx_all_tissues_expression_median',
           'bestWTMatchScore_I', 'bestWTPeptideCount_I', 'rnaseq_gene_expression_quartile',
           'mut_Rank_BA_0', 'mut_Rank_EL_0', 'next_best_BA_mut_ranks', 'next_best_EL_mut_ranks', 'mut_nr_weak_binding_alleles_0',
           'mut_nr_weak_binders_0', 'mut_Rank_Stab_0', 'mut_Rank_Stab_1', 'mut_netchop_Ct_score_0']

cl_vars_16 = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx', 'MIN_MUT_RANK_CI_PRIME',
           'CSCAPE_score', 'GTEx_all_tissues_expression_median',
           'bestWTMatchScore_I', 'bestWTPeptideCount_I', 'rnaseq_gene_expression_quartile',
           'mut_Rank_BA_0', 'mut_Rank_EL_0', 'mut_Rank_BA_1', 'mut_Rank_EL_1', 'mut_nr_weak_binding_alleles_0',
           'mut_nr_weak_binders_0', 'mut_Rank_Stab_0', 'mut_Rank_Stab_1', 'mut_netchop_Ct_score_0']

cl_vars_17 = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx', 'MIN_MUT_RANK_CI_PRIME',
           'CSCAPE_score', 'GTEx_all_tissues_expression_median',
           'bestWTMatchScore_I', 'bestWTPeptideCount_I', 'rnaseq_gene_expression_quartile',
           'mut_Rank_BA_0', 'mut_Rank_EL_0', 'next_best_BA_mut_ranks', 'next_best_EL_mut_ranks', 'mut_nr_weak_binding_alleles_0',
           'mut_nr_weak_binders_0', 'mut_Rank_Stab_0', 'mut_Rank_Stab_1', 'mut_netchop_Ct_score_0', 'mut_netchop_Nt_score_0']

cl_vars_18 = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx', 'mut_aa_coeff_0', 'wt_aa_coeff_0',
           'CSCAPE_score', 'GTEx_all_tissues_expression_median', 'mut_is_binding_pos_0',
           'bestWTMatchScore_I', 'bestWTPeptideCount_I', 'rnaseq_gene_expression_quartile',
           'mut_Rank_BA_0', 'mut_Rank_EL_0', 'next_best_BA_mut_ranks', 'next_best_EL_mut_ranks', 'mut_nr_weak_binding_alleles_0',
           'mut_nr_weak_binders_0', 'mut_Rank_Stab_0', 'mut_Rank_Stab_1', 'mut_netchop_Ct_score_0', 'mut_netchop_Nt_score_0']

cl_vars_19 = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx', 'mut_aa_coeff_0', 'wt_aa_coeff_0',
           'CSCAPE_score', 'GTEx_all_tissues_expression_median',
           'bestWTMatchScore_I', 'bestWTPeptideCount_I', 'rnaseq_gene_expression_quartile',
           'mut_Rank_BA_0', 'mut_Rank_EL_0', 'next_best_BA_mut_ranks', 'next_best_EL_mut_ranks', 'mut_nr_weak_binding_alleles_0',
           'mut_nr_weak_binders_0', 'mut_Rank_Stab_0', 'mut_Rank_Stab_1', 'mut_netchop_Ct_score_0', 'mut_netchop_Nt_score_0']

cl_vars_20 = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx', 'MIN_MUT_RANK_CI_PRIME',
           'CSCAPE_score', 'GTEx_all_tissues_expression_median', 'mut_is_binding_pos_0', 'DAI',
           'bestWTMatchScore_I', 'bestWTPeptideCount_I', 'rnaseq_gene_expression_quartile',
           'mut_Rank_BA_0', 'mut_Rank_EL_0', 'next_best_BA_mut_ranks', 'next_best_EL_mut_ranks', 'mut_nr_weak_binding_alleles_0',
           'mut_nr_weak_binders_0', 'mut_Rank_Stab_0', 'mut_Rank_Stab_1', 'mut_netchop_Ct_score_0', 'mut_netchop_Nt_score_0']

cl_vars = cl_vars_20

print("Features: "+str(cl_vars))

date_time = datetime.now().strftime("%m.%d.%Y-%H.%M.%S")
classifier_file = \
    '/home/localadmin/Priorization/Classifiers/SVM_netMHC_'+date_time+'.sav'
while os.path.isfile(classifier_file):
    date_time = datetime.now().strftime("%m.%d.%Y-%H.%M.%S")
    classifier_file = \
        '/home/localadmin/Priorization/Classifiers/SVM_netMHC_'+date_time+'.sav'

rank_file = '/home/localadmin/Priorization/results/SVM_RBF_rank_data.txt'
f = open(rank_file, "w")
f.close()

print('Rosenberg data class counter: ' + str(Counter(y_R)))

if normalizer is not None:
    X_R = normalizer.fit_transform(df_R[cl_vars])
    X_T1 = normalizer.fit_transform(df_T1[cl_vars])
    X_T2 = normalizer.fit_transform(df_T2[cl_vars])
    X_T3 = normalizer.fit_transform(df_T3[cl_vars])
    X_T4 = normalizer.fit_transform(df_T4[cl_vars])
    X_T8 = normalizer.fit_transform(df_T8[cl_vars])
    X_T9 = normalizer.fit_transform(df_T9[cl_vars])
    X_T12 = normalizer.fit_transform(df_T12[cl_vars])
    X_T16 = normalizer.fit_transform(df_T16[cl_vars])
    X_0YM1 = normalizer.fit_transform(df_0YM1[cl_vars])
    X_058C = normalizer.fit_transform(df_058C[cl_vars])
else:
    X_R = df_R.loc[:, cl_vars].to_numpy()
    X_T1 = df_T1.loc[:, cl_vars].to_numpy()
    X_T2 = df_T2.loc[:, cl_vars].to_numpy()
    X_T3 = df_T3.loc[:, cl_vars].to_numpy()
    X_T4 = df_T4.loc[:, cl_vars].to_numpy()
    X_T8 = df_T8.loc[:, cl_vars].to_numpy()
    X_T9 = df_T9.loc[:, cl_vars].to_numpy()
    X_T12 = df_T12.loc[:, cl_vars].to_numpy()
    X_T16 = df_T16.loc[:, cl_vars].to_numpy()
    X_0YM1 = df_0YM1.loc[:, cl_vars].to_numpy()
    X_058C = df_058C.loc[:, cl_vars].to_numpy()


nr_correct_tot = 0
nr_immuno_tot = 0
nr_missed_tot = 0
nr_correct_adj_tot = 0
nr_immuno_adj_tot = 0
nr_missed_adj_tot = 0

cluster_prob_all = np.array([])
cluster_cnt_all = np.array([])
y_pred_all = np.array([])
score_pred_all = np.array([])
y_all = np.array([])
X_all = None
patient_all = np.array([])
peptide_id_all = np.array([])

best_clfs, best_scores, best_std = \
    optimize_function(X_R, y_R, top100_scorer, create_optim_plots, nr_clfs=nr_clfs, verbose=1,
                      nr_iter=nr_iter)

estimators = []
weights = []
feature_importance = np.zeros(len(cl_vars))
for i in range(len(best_clfs)):
    entry = ('clf_%d' % i, best_clfs[i])
    estimators.append(entry)
    w = best_scores[i] - best_std[i]
    weights.append(w)
    best_clfs[i].fit(X_R, y_R)
    if hasattr(best_clfs[i], 'feature_importances_'):
        feature_importance = feature_importance + np.array(best_clfs[i].feature_importances_)
    elif hasattr(best_clfs[i], 'coef_'):
        feature_importance = feature_importance + np.abs(np.array(best_clfs[i].coef_[0]))

voting_clf = VotingClassifier(estimators=estimators, voting='soft', weights=None)

cluster_probs = UmapClusterProb(X_R,y_R,eps=5.0,verbose=1)

splits = train_test_split(patients, y_R, 0.2)
cv_cnt = 0
cv_clfs = []
cv_weights = []
for train_idx, test_idx in splits:

    cv_cnt += 1

    X_train = X_R[train_idx, :]
    X_test = X_R[test_idx, :]

    y_train = np.array(df_R.loc[train_idx, 'response'])
    y_test = np.array(df_R.loc[test_idx, 'response'])
    patient_test = np.array(df_R.loc[test_idx, 'patient'])
    peptide_id_test = np.array(df_R.loc[test_idx, 'peptide_id'])

    if np.sum(y_test == 1) == 0:
        continue

    voting_clf.fit(X_train, y_train)

    y_all = np.append(y_all, y_test)
    if X_all is None:
        X_all = X_test
    else:
        X_all = np.concatenate((X_all, X_test))

    cluster_prob_test, cluster_cnt_test = cluster_probs.predict_prob(X_test)
    cluster_test = cluster_probs.get_cluster_id(X_test)
    y_pred_all = np.append(y_pred_all, voting_clf.predict_proba(X_test)[:, 1])
    cluster_prob_all = np.append(cluster_prob_all, cluster_prob_test)
    cluster_cnt_all = np.append(cluster_cnt_all, cluster_cnt_test)

    dec_fct = np.full(len(y_test), 0.0)
    for i in range(len(best_clfs)):
        if hasattr(best_clfs[i], 'decision_function'):
            dec_fct += best_clfs[i].decision_function(X_test)
        else:
            dec_fct += best_clfs[i].predict_proba(X_test)[:, 1]

    score_pred_all = np.append(score_pred_all, dec_fct)

    patient_all = np.append(patient_all, patient_test)
    peptide_id_all = np.append(peptide_id_all, peptide_id_test)

    y_pred, nr_correct, nr_immuno, nr_missed, y_pred_adj, nr_correct_adj, nr_immuno_adj, nr_missed_adj = \
        test_classifier_on_patients(voting_clf, patient_test, X_test, y_test, cl_vars, cluster_prob_test,
                                    cluster_cnt_all, cluster_test, 20, verbose=verbose, prob=not use_decision_score,
                                    rank_file=rank_file)

    nr_correct_tot += nr_correct
    nr_immuno_tot += nr_immuno
    nr_missed_tot += nr_missed
    nr_correct_adj_tot += nr_correct_adj
    nr_immuno_adj_tot += nr_immuno_adj
    nr_missed_adj_tot += nr_missed_adj

    print('True positive cnt: %d, immunogenic cnt: %d, missed patient cnt: %d' %(nr_correct, nr_immuno, nr_missed))
    print('Adjusted True positive cnt: %d, immunogenic cnt: %d, missed patient cnt: %d' %
          (nr_correct_adj, nr_immuno_adj, nr_missed_adj))


voting_clf.fit(X_R, y_R)

print("Classifier saved to: " + classifier_file)
pickle.dump(voting_clf, open(classifier_file, 'wb'))
voting_clf = pickle.load(open(classifier_file, 'rb'))

print('True positive cnt: %d, immunogenic cnt: %d, missed patient cnt: %d' %
      (nr_correct_tot, nr_immuno_tot, nr_missed_tot))

print('Adjusted true positive cnt: %d, immunogenic cnt: %d, missed patient cnt: %d' %
      (nr_correct_adj_tot, nr_immuno_adj_tot, nr_missed_adj_tot))

if create_Rosenberg_plots:

    df_test = pd.DataFrame(data=X_all, columns=cl_vars)
    df_test.insert(len(cl_vars), "score", y_pred_all)
    df_test.insert(len(cl_vars) + 1, "response", y_all)
    df_test.insert(len(cl_vars) + 2, "patient", patient_all)
    df_test.insert(len(cl_vars) + 3, "peptide_id", peptide_id_all)
    df_test.insert(len(cl_vars) + 4, "probability", score_pred_all)

    res_file = os.path.join(parameters.get_result_dir(), "Rosenberg_"+class_out_file_tag+".txt")
    df_test.to_csv(res_file, sep = "\t", header=True, index=False)

    df = pd.DataFrame(data=X_all, columns=cl_vars)
    pca_browser = PCAClassifyPeptideBrowser(df, y_all, y_pred_all, patient_all, peptide_id_all, "Rosenberg", True)
    umap_browser = UmapClassifyPeptideBrowser(umaper=cluster_probs, df=df, y=y_all, y_pred=y_pred_all,
                                              patients=patient_all, peptides=peptide_id_all, name="Rosenberg")

    if sum(feature_importance) > 0:
        fig, ax = plt.subplots(figsize=(10, 10))

        x = np.arange(len(feature_importance))
        ax.bar(x, feature_importance)

        ax.set_xticks(x)
        ax.set_xticklabels(cl_vars)
        # plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
        #          rotation_mode="anchor")
        # ax.xticks(x, cl_vars)
        ax.xaxis.set_tick_params(labelsize=4, labelrotation=15)

    plt.show(block=False)


p_T1, cnt_T1 = cluster_probs.predict_prob(X_T1)
cl_T1 = cluster_probs.get_cluster_id(X_T1)
y_pred_T1, nr_correct, nr_immuno, nr_missed, y_pred_adj, nr_correct_adj_tot, nr_immuno_adj_tot, nr_missed_adj_tot = \
    test_classifier_on_patients(voting_clf, np.full(df_T1.shape[0],'TESLA1'), X_T1, y_T1, cl_vars, p_T1, cnt_T1, cl_T1,
                                20, verbose=verbose, prob=True)

p_T2, cnt_T2 = cluster_probs.predict_prob(X_T2)
cl_T2 = cluster_probs.get_cluster_id(X_T2)
y_pred_T2, nr_correct, nr_immuno, nr_missed, y_pred_adj, nr_correct_adj_tot, nr_immuno_adj_tot, nr_missed_adj_tot = \
    test_classifier_on_patients(voting_clf, np.full(df_T2.shape[0], 'TESLA2'), X_T2, y_T2, cl_vars, p_T2, cnt_T2, cl_T2,
                                20, verbose=verbose, prob=True)

p_T3, cnt_T3 = cluster_probs.predict_prob(X_T3)
cl_T3 = cluster_probs.get_cluster_id(X_T3)
y_pred_T3, nr_correct, nr_immuno, nr_missed, y_pred_adj, nr_correct_adj_tot, nr_immuno_adj_tot, nr_missed_adj_tot = \
    test_classifier_on_patients(voting_clf, np.full(df_T3.shape[0], 'TESLA3'), X_T3, y_T3, cl_vars, p_T3, cnt_T3, cl_T3,
                                20, verbose=verbose, prob=True)

p_T4, cnt_T4 = cluster_probs.predict_prob(X_T4)
cl_T4 = cluster_probs.get_cluster_id(X_T4)
y_pred_T4, nr_correct, nr_immuno, nr_missed, y_pred_adj, nr_correct_adj_tot, nr_immuno_adj_tot, nr_missed_adj_tot = \
    test_classifier_on_patients(voting_clf, np.full(df_T4.shape[0], 'TESLA4'), X_T4, y_T4, cl_vars, p_T4, cnt_T4, cl_T4,
                                20, verbose=verbose, prob=True)

p_T8, cnt_T8 = cluster_probs.predict_prob(X_T8)
cl_T8 = cluster_probs.get_cluster_id(X_T8)
y_pred_T8, nr_correct, nr_immuno, nr_missed, y_pred_adj, nr_correct_adj_tot, nr_immuno_adj_tot, nr_missed_adj_tot = \
    test_classifier_on_patients(voting_clf, np.full(df_T8.shape[0], 'TESLA8'), X_T8, y_T8, cl_vars, p_T8, cnt_T8, cl_T8,
                                20, verbose=verbose, prob=True)

p_T9, cnt_T9 = cluster_probs.predict_prob(X_T9)
cl_T9 = cluster_probs.get_cluster_id(X_T9)
y_pred_T9, nr_correct, nr_immuno, nr_missed, y_pred_adj, nr_correct_adj_tot, nr_immuno_adj_tot, nr_missed_adj_tot = \
    test_classifier_on_patients(voting_clf, np.full(df_T9.shape[0], 'TESLA9'), X_T9, y_T9, cl_vars, p_T9, cnt_T9, cl_T9,
                                20, verbose=verbose, prob=True)

p_T12, cnt_T12 = cluster_probs.predict_prob(X_T12)
cl_T12 = cluster_probs.get_cluster_id(X_T12)
y_pred_T12, nr_correct, nr_immuno, nr_missed, y_pred_adj, nr_correct_adj_tot, nr_immuno_adj_tot, nr_missed_adj_tot = \
    test_classifier_on_patients(voting_clf, np.full(df_T12.shape[0], 'TESLA12'), X_T12, y_T12, cl_vars, p_T12, cnt_T12, cl_T12,
                                20, verbose=verbose, prob=True)

p_T16, cnt_T16 = cluster_probs.predict_prob(X_T16)
cl_T16 = cluster_probs.get_cluster_id(X_T16)
y_pred_T16, nr_correct, nr_immuno, nr_missed, y_pred_adj, nr_correct_adj_tot, nr_immuno_adj_tot, nr_missed_adj_tot = \
    test_classifier_on_patients(voting_clf, np.full(df_T16.shape[0], 'TESLA16'), X_T16, y_T16, cl_vars, p_T16, cnt_T16, cl_T16,
                                20, verbose=verbose, prob=True)


p_0YM1, cnt_0MY1 = cluster_probs.predict_prob(X_0YM1)
cl_0YM1 = cluster_probs.get_cluster_id(X_0YM1)
y_pred_0YM1, nr_correct, nr_immuno, nr_missed, y_pred_adj, nr_correct_adj_tot, nr_immuno_adj_tot, nr_missed_adj_tot = \
    test_classifier_on_patients(voting_clf, np.full(df_0YM1.shape[0], '0YM1'), X_0YM1, y_0YM1, cl_vars, p_0YM1, cnt_0MY1,
                                cl_0YM1, 20, verbose=verbose, prob=True)

p_058C, cnt_058C = cluster_probs.predict_prob(X_058C)
cl_058C = cluster_probs.get_cluster_id(X_058C)
y_pred_058C, nr_correct, nr_immuno, nr_missed, y_pred_adj, nr_correct_adj_tot, nr_immuno_adj_tot, nr_missed_adj_tot = \
    test_classifier_on_patients(voting_clf, np.full(df_058C.shape[0], '058C'), X_058C, y_058C, cl_vars, p_058C,
                                cnt_058C, cl_058C, 20, verbose=verbose, prob=True)


if create_Tesla_plots:
    create_report(df_T1, X_T1, cl_vars, y_pred_T1, y_T1, "TESLA1", class_out_file_tag, cluster_probs)
    create_report(df_T2, X_T2, cl_vars, y_pred_T2, y_T2, "TESLA2", class_out_file_tag, cluster_probs)
    create_report(df_T3, X_T3, cl_vars, y_pred_T3, y_T3, "TESLA3", class_out_file_tag, cluster_probs)
    create_report(df_T4, X_T4, cl_vars, y_pred_T4, y_T4, "TESLA4", class_out_file_tag, cluster_probs)
    create_report(df_T8, X_T8, cl_vars, y_pred_T8, y_T8, "TESLA8", class_out_file_tag, cluster_probs)
    create_report(df_T9, X_T9, cl_vars, y_pred_T9, y_T9, "TESLA9", class_out_file_tag, cluster_probs)
    create_report(df_T12, X_T12, cl_vars, y_pred_T12, y_T12, "TESLA12", class_out_file_tag, cluster_probs)
    create_report(df_T16, X_T16, cl_vars, y_pred_T16, y_T16, "TESLA16", class_out_file_tag, cluster_probs)


if create_0YM1_plots:
    create_report(df_0YM1, X_0YM1, cl_vars, y_pred_0YM1, y_0YM1, "0YM1", class_out_file_tag, cluster_probs)

if create_058C_plots:
    create_report(df_058C, X_058C, cl_vars, y_pred_058C, y_058C, "058C", class_out_file_tag, cluster_probs)


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



