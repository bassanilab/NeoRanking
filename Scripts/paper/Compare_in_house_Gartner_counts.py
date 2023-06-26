from DataWrangling.DataTransformer import DataTransformer
from DataWrangling.RosenbergImmunogenicityAnnotatorLong import *
from DataWrangling.RosenbergImmunogenicityAnnotatorShort import *
from Utils.DataManager import *


def contains_mutation(inhouse_data, mutant_sequence):
    return any(inhouse_data['mutant_seq'] == mutant_sequence)


def get_overlapping_mutation_count(gartner_data, inhouse_data):
    return sum(gartner_data.apply(lambda row: contains_mutation(inhouse_data, row['Mut Epitope']), axis=1))


rosenberg_info_train_file, rosenberg_info_test_file = GlobalParameters().get_gartner_info_files(peptide_type='long')
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

inhouse_info = []
for p in patients:
    if p in mgr.get_valid_patients(peptide_type='long') and p in \
            set.union(annotator_long.get_patients('nci_test'), annotator_long.get_patients('nci_train')):
        inhouse_data = annotator_long.annotate_patient(p)
        inhouse_data = inhouse_data.astype({'mutant_seq': 'string'})

        gartner_data_sel = gartner_data.loc[gartner_data['ID'] == p, ]

        gartner_data_sel_SNV = gartner_data_sel.loc[gartner_data_sel['Mutation type'] == 'nonsynonymous SNV', :]
        inhouse_data_SNV = inhouse_data.loc[inhouse_data['mutation_type'] == 'SNV', :]

        d = {'Patient': p,
             'NCI_train_mut_seq_all in-house': inhouse_data.shape[0],
             'NCI_train_mut_seq_tested in-house': sum(np.logical_or(inhouse_data['response_type'] == 'negative',
                                                     inhouse_data.apply(lambda r: 'CD8' in r['response_type'], axis=1))),
             'NCI_train_mut_seq_imm in-house SNV':
                 sum(np.logical_and(inhouse_data['mutation_type'] == 'SNV',
                                    inhouse_data.apply(lambda r: 'CD8' in r['response_type'], axis=1))),
             'NCI_train_mut_seq_imm in-house': sum(inhouse_data.apply(lambda r: 'CD8' in r['response_type'], axis=1)),
             'NCI_train_mut_seq in-house SNV': sum(inhouse_data['mutation_type'] == 'SNV'),
             'NCI_train_mut_seq in-house InDel': sum(np.logical_or(inhouse_data['mutation_type'] == 'INSERTION',
                                                    inhouse_data['mutation_type'] == 'DELETION')),
             'NCI_train_mut_seq in-house FS': sum(inhouse_data['mutation_type'] == 'FSS'),
             'NCI_train_mut_seq_all Gartner': gartner_data_sel.shape[0],
             'NCI_train_mut_seq_tested Gartner': sum(np.logical_or(inhouse_data['response_type'] == 'negative',
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
             'NCI_train_mut_seq_all in-house/Gartner overlap': get_overlapping_mutation_count(gartner_data_sel, inhouse_data),
             'NCI_train_mut_seq in-house/Gartner SNV overlap': get_overlapping_mutation_count(gartner_data_sel_SNV, inhouse_data_SNV)}

        result_df = result_df.append(pd.Series(d), ignore_index=True)

        result_df['NCI_train_mut_seq_all in-house/Gartner overlap ratio'] = \
            result_df['NCI_train_mut_seq_all in-house/Gartner overlap']/result_df['NCI_train_mut_seq_all Gartner']

        result_df['NCI_train_mut_seq in-house/Gartner SNV overlap ratio'] = \
            result_df['NCI_train_mut_seq in-house/Gartner SNV overlap']/result_df['NCI_train_mut_seq Gartner SNV']

        order = ['Patient', 'NCI_train_mut_seq_all in-house', 'NCI_train_mut_seq_all Gartner',
                 'NCI_train_mut_seq_all in-house/Gartner overlap', 'NCI_train_mut_seq_all in-house/Gartner overlap ratio',
                 'NCI_train_mut_seq_tested in-house', 'NCI_train_mut_seq_tested Gartner',
                 'NCI_train_mut_seq_imm in-house', 'NCI_train_mut_seq_imm Gartner',
                 'NCI_train_mut_seq_imm in-house SNV', 'NCI_train_mut_seq_imm Gartner SNV',
                 'NCI_train_mut_seq in-house SNV', 'NCI_train_mut_seq Gartner SNV',
                 'NCI_train_mut_seq in-house/Gartner SNV overlap', 'NCI_train_mut_seq in-house/Gartner SNV overlap ratio',
                 'NCI_train_mut_seq in-house InDel', 'NCI_train_mut_seq Gartner InDel',
                 'NCI_train_mut_seq in-house FS', 'NCI_train_mut_seq Gartner FS']
        result_df = result_df[order]

result_df.to_csv(os.path.join(GlobalParameters().get_plot_dir(), 'Compare_in-house_Gartner_Counts.txt'),
                 header=True, index=False, sep="\t")
