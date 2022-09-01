from DataWrangling.DataLoader import DataLoader
from DataWrangling.RosenbergImmunogenicityAnnotatorLong import *
from DataWrangling.RosenbergImmunogenicityAnnotatorShort import *
from Utils.DataManager import *


def contains_mutation(neodisc_data, mutant_sequence):
    return any(neodisc_data['mutant_seq'] == mutant_sequence)


def get_overlapping_mutation_count(gartner_data, neodisc_data):
    return sum(gartner_data.apply(lambda row: contains_mutation(neodisc_data, row['Mut Epitope']), axis=1))


rosenberg_info_train_file, rosenberg_info_test_file = Parameters().get_gartner_info_files(peptide_type='long')
data_train = pd.read_csv(rosenberg_info_train_file, sep="\t", header=0)
data_test = pd.read_csv(rosenberg_info_test_file, sep="\t", header=0)

data_train.loc[data_train['Screening Status'] == '-', 'Screening Status'] = 'negative'
data_train.loc[data_train['Screening Status'] == 'unscreened', 'Screening Status'] = 'not_tested'

data_test.loc[data_test['Screening Status'] == '0', 'Screening Status'] = 'negative'
data_test.loc[data_test['Screening Status'] == '1', 'Screening Status'] = 'CD8'
data_test.loc[data_test['Screening Status'] == 'unscreened', 'Screening Status'] = 'not_tested'

gartner_data = pd.concat([data_test, data_train], ignore_index=True, join='inner')
gartner_data = gartner_data.astype({'ID': str, 'Mut Epitope': str})

patients = set(gartner_data['ID'].unique())

result_df = pd.DataFrame()
mgr = DataManager(rna_seq=True, netmhc=True)
annotator_long = RosenbergImmunogenicityAnnotatorLong(mgr)

neodisc_info = []
for p in patients:
    if p in mgr.get_valid_patients(peptide_type='long') and p in \
            set.union(annotator_long.get_patients('gartner_test'), annotator_long.get_patients('gartner_train')):
        neodisc_data = annotator_long.annotate_response_types(p)
        neodisc_data = neodisc_data.astype({'mutant_seq': 'string'})

        gartner_data_sel = gartner_data.loc[gartner_data['ID'] == p,]

        gartner_data_sel_SNV = gartner_data_sel.loc[gartner_data_sel['Mutation type'] == 'nonsynonymous SNV', :]
        neodisc_data_SNV = neodisc_data.loc[neodisc_data['mutation_type'] == 'SNV', :]

        d = {'Patient': p,
             'NeoDisc mut_all': neodisc_data.shape[0],
             'NeoDisc mut_tested': sum(np.logical_or(neodisc_data['response_type'] == 'negative',
                                                     neodisc_data['response_type'] == 'CD8')),
             'NeoDisc mut_imm': sum(neodisc_data['response_type'] == 'CD8'),
             'NeoDisc SNV mut': sum(neodisc_data['mutation_type'] == 'SNV'),
             'NeoDisc InDel mut': sum(np.logical_or(neodisc_data['mutation_type'] == 'INSERTION',
                                                    neodisc_data['mutation_type'] == 'DELETION')),
             'NeoDisc FS mut': sum(neodisc_data['mutation_type'] == 'FSS'),
             'Gartner mut_all': gartner_data_sel.shape[0],
             'Gartner mut_tested': sum(np.logical_or(neodisc_data['response_type'] == 'negative',
                                                     neodisc_data['response_type'] == 'CD8')),
             'Gartner mut_imm':
                 sum(gartner_data_sel['Screening Status'] == 'CD8'),
             'Gartner SNV mut': sum(gartner_data_sel['Mutation type'] == 'nonsynonymous SNV'),
             'Gartner InDel mut': sum(np.logical_or(gartner_data_sel['Mutation type'] == 'nonframeshift insertion',
                                                    gartner_data_sel['Mutation type'] == 'nonframeshift deletion')),
             'Gartner FS mut': sum(np.logical_or(gartner_data_sel['Mutation type'] == 'frameshift insertion',
                                                 gartner_data_sel['Mutation type'] == 'frameshift deletion')),
             'NeoDisc-Gartner mut_all overlap': get_overlapping_mutation_count(gartner_data_sel, neodisc_data),
             'NeoDisc-Gartner SNV mut overlap': get_overlapping_mutation_count(gartner_data_sel_SNV, neodisc_data_SNV)}

        result_df = result_df.append(pd.Series(d), ignore_index=True)

        result_df['NeoDisc-Gartner mut_all overlap ratio'] = \
            result_df['NeoDisc-Gartner mut_all overlap']/result_df['Gartner mut_all']

        result_df['NeoDisc-Gartner SNV mut overlap ratio'] = \
            result_df['NeoDisc-Gartner SNV mut overlap']/result_df['Gartner SNV mut']

        order = ['Patient', 'NeoDisc mut_all', 'Gartner mut_all', 'NeoDisc-Gartner mut_all overlap',
                 'NeoDisc-Gartner mut_all overlap ratio',
                 'NeoDisc mut_tested', 'Gartner mut_tested', 'NeoDisc mut_imm', 'Gartner mut_imm',
                 'NeoDisc SNV mut', 'Gartner SNV mut', 'NeoDisc-Gartner SNV mut overlap',
                 'NeoDisc-Gartner SNV mut overlap ratio',
                 'NeoDisc InDel mut', 'Gartner InDel mut', 'NeoDisc FS mut', 'Gartner FS mut']
        result_df = result_df[order]

result_df.to_csv(os.path.join(Parameters().get_plot_dir(), 'Compare_NeoDisc_Gartner_Counts.txt'),
                 header=True, index=False, sep="\t")
