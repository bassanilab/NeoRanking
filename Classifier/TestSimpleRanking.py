import argparse
import re
from catboost import CatBoostClassifier

from DataWrangling.DataTransformer import *
from Classifier.PrioritizationLearner import *
from Utils.Util_fct import *
from sklearn.ensemble import VotingClassifier

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-d', '--result_dir', type=str, default=GlobalParameters().get_pickle_dir(),
                    help='directory for result files')
parser.add_argument('-sr', '--simple_ranking', type=str, default='netmhc',
                    help='Simple ranking method (mixmhc or MuPeXI')
parser.add_argument('-te', '--patients_test', type=str, default='test', help='patient ids for test set')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-r', '--max_rank', type=int, default=20,
                    help='Maximal rank for predicted immunogenic considered correct')
parser.add_argument('-mi', '--min_nr_immuno', type=int, default=1, help='Minimum nr of immunogenic mutations in sample')
parser.add_argument('-a', '--alpha', type=float, default=0.05, help='Coefficient alpha in score function')
parser.add_argument('-cat', '--cat_encoder', type=str, default='float', help='convert categories to numbers')
parser.add_argument('-pt', '--peptide_type', type=str, default='short', help='Peptide type (long or short)')
parser.add_argument('-eg', '--excluded_genes', type=str, nargs='+', help='genes excluded from prioritization')
parser.add_argument('-mrn', '--max_rank_netmhc', type=int, default=20000, help='Maximal netmhc rank of short peptide')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))


patients_test = \
    get_valid_patients(dataset=args.patients_test, peptide_type=args.peptide_type) \
        if args.patients_test and len(args.patients_test) > 0 else get_valid_patients(peptide_type=args.peptide_type)

mgr = DataManager()
patients_test = sorted(patients_test.intersection(mgr.get_immunogenic_patients(args.peptide_type)))

encodings = read_cat_encodings('NCI_train', args.peptide_type)

data_loader_test = DataTransformer(transformer=DataTransformer(),
                                   features=['mutant_rank', 'wt_best_rank', 'rnaseq_TPM', 'rnaseq_ref_support'
                                        'mutant_rank_netMHCpan', 'wt_rank_netMHCpan'],
                                   mutation_types=['SNV'], response_types=['CD8', 'CD4/CD8', 'negative', 'not_tested'],
                                   immunogenic=['CD8', 'CD4/CD8'], min_nr_immuno=1, cat_type='int', max_netmhc_rank=20,
                                   cat_encoders=encodings, excluded_genes=args.excluded_genes)


def sum_rank_correct(ranks, alpha):
    return np.sum(np.exp(np.multiply(-alpha, ranks)))


def write_results(patient, data, report_file, max_rank):
    if report_file and os.path.getsize(report_file.name) == 0:
        report_file.write("Patient\tNr_correct_top{0}\tNr_immunogenic\tMax_rank\tNr_peptides\tClf_score\t"
                          "CD8_ranks\tCD8_peptide_idx\tCD8_mut_seqs\tCD8_genes\n".format(max_rank))

    r = np.where(data['response'] == 1)[0]
    nr_correct = sum(r < max_rank)
    nr_immuno = len(r)
    score = sum_rank_correct(r, args.alpha)
    sort_idx = np.argsort(r)
    ranks_str = ",".join(["{0:.0f}".format(np.floor(r+1)) for r in r[sort_idx]])
    peptide_ids = data.loc[data['response'] == 1, 'peptide_id'].to_numpy()
    peptide_id_str = ",".join(["{0}".format(s) for s in peptide_ids[sort_idx]])
    mut_seqs = data.loc[data['response'] == 1, 'mutant_seq'].to_numpy()
    mut_seqs_str = ",".join(["{0}".format(s) for s in mut_seqs[sort_idx]])
    genes = data.loc[data['response'] == 1, 'gene'].to_numpy()
    gene_str = ",".join(["{0}".format(s) for s in genes[sort_idx]])

    if report_file:
        report_file.write("%s\t%d\t%d\t%d\t%d\t%f\t%s\t%s\t%s\t%s\n" %
                          (patient, nr_correct, nr_immuno, np.min((max_rank, data.shape[0])), data.shape[0], score,
                           ranks_str, peptide_id_str, mut_seqs_str, gene_str))


def MuPeXI_score(row):
    L1 = 1.0/(1 + np.exp(5*(row['mutant_rank']-2)))
    L2 = 1.0/(1 + np.exp(5*(row['wt_best_rank']-2)))
    A = row['rnaseq_ref_support']/100
#    return L1*A*np.tanh(row['rnaseq_TPM']/200)*(1-L2/2)
#    return L1*A*row['rnaseq_TPM']#*(1-L2/2)
#    return L1*row['rnaseq_TPM']#*(1-L2/2)
    return L1*row['rnaseq_TPM']*(1-L2/2)


def sort(method, data):
    if method == 'mixmhc':
        sorted_data = data.copy()
        sorted_data.loc[:, 'rnaseq_TPM'] *= -1
        sorted_data = data.sort_values(by=['mutant_rank', 'rnaseq_TPM'], ascending=True)
        sorted_data.loc[:, 'rnaseq_TPM'] *= -1
        return sorted_data

    if method == 'netmhc':
        sorted_data = data.copy()
        sorted_data.loc[:, 'rnaseq_TPM'] *= -1
        sorted_data = data.sort_values(by=['mutant_rank_netMHCpan', 'rnaseq_TPM'], ascending=True)
        sorted_data.loc[:, 'rnaseq_TPM'] *= -1
        return sorted_data

    if method == 'MuPeXI':
        scores = data.apply(lambda row: MuPeXI_score(row), axis=1)
        data['MuPeXI'] = scores
        sorted_data = data.sort_values(by=['MuPeXI'], ascending=False)
        return sorted_data


def get_learner(classifier_name, x):
    optimizationParams = OptimizationParams(args.alpha, cat_idx=data_loader_test.get_categorical_idx(x),
                                            cat_dims=data_loader_test.get_categorical_dim(),
                                            input_shape=[len(args.features)])

    return PrioritizationLearner(classifier_name, args.scorer, optimizationParams, verbose=args.verbose)


result_file_name = \
    os.path.join(args.result_dir, "SimpleRanking_{0}_{1}_test.txt".format(args.peptide_type, args.simple_ranking))

with open(result_file_name, mode='w') as result_file:

    for p in patients_test:
        data_test, X_test, y_test = data_loader_test.load_patients(p, args.input_file_tag, args.peptide_type, verbose=True)
        if data_test is None:
            continue

        data_test = sort(args.simple_ranking, data_test)
        write_results(p, data_test, result_file, args.max_rank)

