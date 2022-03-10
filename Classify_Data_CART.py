from sklearn.ensemble import VotingClassifier
from sklearn.metrics import make_scorer
from Classify_functions import *
import warnings
from sklearn.exceptions import UndefinedMetricWarning
from sklearn.tree import export_graphviz
import pydotplus
import os

warnings.filterwarnings(action='ignore', category=UndefinedMetricWarning)

create_class_plots = False
create_optim_plots = False
verbose = 0
nr_clfs = 1
nr_iter = 10

top100_scorer = make_scorer(nr_correct_top100, needs_threshold = True)
top100_scorer = make_scorer(sum_prob_correct, needs_threshold = True)

data_file = \
    '~/Documents/data/Prediction/Rosenberg_data/25.08.2020/RosenbergSamplesSummary_Preds_ipMSDB.txt'

res_dir = '/Users/markusmueller/Documents/data/Prediction/Plots_Tesla_Rosenberg'
res_file = os.path.join(res_dir,'RosenbergSamplesSummary_CART_results.txt')

df, y = load_data(data_file)
df.insert(0, "response", y)

print('Class counter: ' + str(Counter(y)))

patients = np.array(df.patient)

cl_vars = ['rnaseq_ref_support','Sample_Tissue_expression_GTEx','MIN_MUT_RANK_CI_MIXMHC','WT_BEST_RANK_CI_MIXMHC',
           'COUNT_MUT_RANK_CI_MIXMHC','WT_RANK_CI_MIXMHC','MIN_MUT_RANK_CI_PRIME','WT_BEST_RANK_CI_PRIME',
           'COUNT_MUT_RANK_CI_PRIME','WT_RANK_CI_PRIME','CSCAPE_score','GTEx_all_tissues_expression_median','bestWTMatchScore_I',
           'bestWTMatchOverlap_I','bestWTPeptideCount_I','bestWT_Cancer_ProteinCnt_I','rnaseq_gene_expression_quartile']

X = df[cl_vars]

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

splits = train_test_split(patients, y, 0.2)
cv_cnt = 0
cv_carts = []
cv_weights = []
for train_idx, test_idx in splits:

    cv_cnt += 1
    train_idx = df.index[train_idx]
    test_idx = df.index[test_idx]

    X_train = df.loc[train_idx,cl_vars].to_numpy()
    X_test = df.loc[test_idx,cl_vars].to_numpy()

    #X_org_test = df.loc[test_idx,cl_vars]
#    X_org_test = df.loc[test_idx,]
    y_train = np.array(df.loc[train_idx,'response'])
    y_test = np.array(df.loc[test_idx,'response'])
    p_test = np.array(df.loc[test_idx,'patient'])
    id_test = np.array(df.loc[test_idx,'peptide_id'])

    if np.sum(y_test == 1) == 0:
        continue

    best_cart, best_scores, best_std = \
        optimize_cart_rdn(X_train, y_train, top100_scorer, create_optim_plots, nr_clfs=nr_clfs, verbose=0,
                          nr_iter=nr_iter)

    estimators = []
    weights = []
    for i in range(len(best_cart)):
        entry = ('svm_%d_%d' %(cv_cnt, i), best_cart[i])
        estimators.append(entry)
        cv_carts.append(entry)
        w = best_scores[i]-best_std[i]
        weights.append(w)
        cv_weights.append(w)
        best_cart[i].fit(X_train, y_train)

    voting_cart = VotingClassifier(estimators = estimators, voting='soft',weights=None)
    voting_cart.fit(X_train, y_train)

    y_all = np.append(y_all,y_test)
    if X_all is None:
        X_all = X_test
    else:
        X_all = np.concatenate((X_all, X_test))

    y_pred_all = np.append(y_pred_all, voting_cart.predict_proba(X_test)[:,1])
    score_pred_all = y_pred_all
    cv_layer = np.append(cv_layer, np.full(len(p_test),cv_cnt))
    p_all = np.append(p_all, p_test)
    id_all = np.append(id_all, id_test)

    y_pred, nr_correct, nr_immuno, nr_missed = \
        test_classifier_on_patients(voting_cart, p_test, X_test, y_test, 20, verbose=1, prob=True)

    nr_correct_tot += nr_correct
    nr_immuno_tot += nr_immuno
    nr_missed_tot +=  nr_missed

    dot_data = export_graphviz(best_cart[0], out_file=None, class_names=['negative',"CD8"],
                               feature_names=cl_vars,rounded=True, filled=True)
    graph = pydotplus.graphviz.graph_from_dot_data(dot_data)
    out_file = os.path.join(res_dir,'CART_decision_tree_'+str(cv_cnt)+'.png')
    graph.write_png(out_file)

    print('True positive cnt: %d, immunogenic cnt: %d, missed patient cnt: %d' %(nr_correct, nr_immuno, nr_missed))


print('True positive cnt: %d, immunogenic cnt: %d, missed patient cnt: %d' % (nr_correct_tot, nr_immuno_tot,
                                                                              nr_missed_tot))
voting_cart = VotingClassifier(estimators=cv_carts, voting='soft', weights=cv_weights)
voting_cart.fit(df[cl_vars], y)

data_file = \
    '~/Documents/data/Prediction/Tesla/PredictionFiles/26.08.2020/TESLA1.long.rnaseq.sample.tumorcontent.immunogenic.ipMSDB.txt'

df_T1, y_T1 = load_data(data_file)

nr_correct, nr_immuno, nr_missed = test_classifier_on_patients(voting_cart, np.full(df_T1.shape[0],'TESLA1'), df_T1[cl_vars], y_T1, 20,
                                                               verbose=1, prob=True)

data_file = \
    '~/Documents/data/Prediction/Tesla/PredictionFiles/26.08.2020/TESLA3.long.rnaseq.sample.tumorcontent.immunogenic.ipMSDB.txt'

df_T3, y_T3 = load_data(data_file)

nr_correct, nr_immuno, nr_missed = test_classifier_on_patients(voting_cart, np.full(df_T3.shape[0],'TESLA3'), df_T3[cl_vars], y_T3, 20,
                                                               verbose=1, prob=True)

if create_class_plots:

    df_test = pd.DataFrame(data=X_all, columns=cl_vars)
    df_test.insert(len(cl_vars), "score", y_pred_all)
    df_test.insert(len(cl_vars) + 1, "response", y_all)
    df_test.insert(len(cl_vars) + 2, "patient", p_all)
    df_test.insert(len(cl_vars) + 3, "peptide_id", id_all)
    df_test.insert(len(cl_vars) + 4, "probability", score_pred_all)

    df_test.to_csv(res_file, sep = "\t", header=True, index=False)

    plot_class_results(df_test)
    plt.show()


# gamma : if overfitting -> decrease gamma. kernel(x','s) = exp[-gamma*(x-s)^2] : the higher gamma the more local the kernel
# C :if overfitting -> decrease C. min 1/2*w*w^t + C*sum_i(eta_i) : eta_i slack variable = (0 if class true',' dist to margin otherwise)
# class_weight : eg. class_weight={0: 1',' 1: 2}. C is replaced by C*class_weight[i]

# best so far {'C': 1.5',' 'class_weight': {1: 9}',' 'gamma': 1}
# Best params: {'C': 0.2750980053342904',' 'class_weight': {1: 11}',' 'gamma': 2.106785011385351',' 'kernel': 'rbf'}',' best precision: 0.11773887232449073
# Best params: {'C': 0.23091500814262222',' 'class_weight': {1: 10}',' 'gamma': 0.66825290705259',' 'kernel': 'rbf'}',' best precision: 0.2000549224220788
# 2885  C:0.213775       Gamma:0.73727      cl_weights:{1: 10}                      mean_precision:0.200055           0.400041                   1
# Best params: {'C': 0.06175348288740734',' 'class_weight': {1: 10}',' 'gamma': 3.609993861334124',' 'kernel': 'rbf'}',' nr_correct_top20: 2.0



