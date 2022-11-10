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
            set.union(annotator_long.get_patients('nci_test'), annotator_long.get_patients('nci_train')):
        neodisc_data = annotator_long.annotate_patient(p)
        neodisc_data = neodisc_data.astype({'mutant_seq': 'string'})

        gartner_data_sel = gartner_data.loc[gartner_data['ID'] == p, ]

        gartner_data_sel_SNV = gartner_data_sel.loc[gartner_data_sel['Mutation type'] == 'nonsynonymous SNV', :]
        neodisc_data_SNV = neodisc_data.loc[neodisc_data['mutation_type'] == 'SNV', :]

        d = {'Patient': p,
             'NCI_train_mut_seq_all NeoDisc': neodisc_data.shape[0],
             'NCI_train_mut_seq_tested NeoDisc': sum(np.logical_or(neodisc_data['response_type'] == 'negative',
                                                     neodisc_data.apply(lambda r: 'CD8' in r['response_type'], axis=1))),
             'NCI_train_mut_seq_imm NeoDisc SNV':
                 sum(np.logical_and(neodisc_data['mutation_type'] == 'SNV',
                                    neodisc_data.apply(lambda r: 'CD8' in r['response_type'], axis=1))),
             'NCI_train_mut_seq_imm NeoDisc': sum(neodisc_data.apply(lambda r: 'CD8' in r['response_type'], axis=1)),
             'NCI_train_mut_seq NeoDisc SNV': sum(neodisc_data['mutation_type'] == 'SNV'),
             'NCI_train_mut_seq NeoDisc InDel': sum(np.logical_or(neodisc_data['mutation_type'] == 'INSERTION',
                                                    neodisc_data['mutation_type'] == 'DELETION')),
             'NCI_train_mut_seq NeoDisc FS': sum(neodisc_data['mutation_type'] == 'FSS'),
             'NCI_train_mut_seq_all Gartner': gartner_data_sel.shape[0],
             'NCI_train_mut_seq_tested Gartner': sum(np.logical_or(neodisc_data['response_type'] == 'negative',
                                                     gartner_data_sel['Screening Status'] == 'CD8')),
             'NCI_train_mut_seq_imm Gartner': sum(gartner_data_sel['Screening Status'] == 'CD8'),
             'NCI_train_mut_seq_imm Gartner SNV':
                 sum(np.logical_and(gartner_data_sel['Mutation type'] == 'nonsynonymous SNV',
                                    gartner_data_sel['Screening Status'] == 'CD8')),
             'NCI_train_mut_seq Gartner SNV': sum(gartner_data_sel['Mutation type'] == 'nonsynonymous SNV'),
             'NCI_train_mut_seq Gartner InDel': sum(np.logical_or(gartner_data_sel['Mutation type'] == 'nonframeshift insertion',
                                                    gartner_data_sel['Mutation type'] == 'nonframeshift deletion')),
             'NCI_train_mut_seq Gartner FS': sum(np.logical_or(gartner_data_sel['Mutation type'] == 'frameshift insertion',
                                                 gartner_data_sel['Mutation type'] == 'frameshift deletion')),
             'NCI_train_mut_seq_all NeoDisc-Gartner overlap': get_overlapping_mutation_count(gartner_data_sel, neodisc_data),
             'NCI_train_mut_seq NeoDisc-Gartner SNV overlap': get_overlapping_mutation_count(gartner_data_sel_SNV, neodisc_data_SNV)}

        result_df = result_df.append(pd.Series(d), ignore_index=True)

        result_df['NCI_train_mut_seq_all NeoDisc-Gartner overlap ratio'] = \
            result_df['NCI_train_mut_seq_all NeoDisc-Gartner overlap']/result_df['NCI_train_mut_seq_all Gartner']

        result_df['NCI_train_mut_seq NeoDisc-Gartner SNV overlap ratio'] = \
            result_df['NCI_train_mut_seq NeoDisc-Gartner SNV overlap']/result_df['NCI_train_mut_seq Gartner SNV']

        order = ['Patient', 'NCI_train_mut_seq_all NeoDisc', 'NCI_train_mut_seq_all Gartner',
                 'NCI_train_mut_seq_all NeoDisc-Gartner overlap', 'NCI_train_mut_seq_all NeoDisc-Gartner overlap ratio',
                 'NCI_train_mut_seq_tested NeoDisc', 'NCI_train_mut_seq_tested Gartner',
                 'NCI_train_mut_seq_imm NeoDisc', 'NCI_train_mut_seq_imm Gartner',
                 'NCI_train_mut_seq_imm NeoDisc SNV', 'NCI_train_mut_seq_imm Gartner SNV',
                 'NCI_train_mut_seq NeoDisc SNV', 'NCI_train_mut_seq Gartner SNV',
                 'NCI_train_mut_seq NeoDisc-Gartner SNV overlap', 'NCI_train_mut_seq NeoDisc-Gartner SNV overlap ratio',
                 'NCI_train_mut_seq NeoDisc InDel', 'NCI_train_mut_seq Gartner InDel',
                 'NCI_train_mut_seq NeoDisc FS', 'NCI_train_mut_seq Gartner FS']
        result_df = result_df[order]

result_df.to_csv(os.path.join(Parameters().get_plot_dir(), 'Compare_NeoDisc_Gartner_Counts.txt'),
                 header=True, index=False, sep="\t")
