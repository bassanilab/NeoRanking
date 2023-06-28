from DataWrangling.DataTransformer import DataTransformer
from Utils.DataManager import *


def contains_mutation(inhouse_data, mutant_sequence):
    return any(inhouse_data['mutant_seq'] == mutant_sequence)


def get_overlapping_mutation_count(gartner_data, inhouse_data):
    return sum(gartner_data.apply(lambda row: contains_mutation(inhouse_data, row['Mut Epitope']), axis=1))


if __name__ == "__main__":

    gartner_info_train_file = GlobalParameters.gartner_nmer_train_file
    gartner_info_test_file = GlobalParameters.gartner_nmer_test_file
    data_train = pd.read_csv(gartner_info_train_file, sep="\t", header=0)
    data_test = pd.read_csv(gartner_info_test_file, sep="\t", header=0)

    data_train.loc[data_train['Screening Status'] == '-', 'Screening Status'] = 'negative'
    data_train.loc[data_train['Screening Status'] == 'unscreened', 'Screening Status'] = 'not_tested'

    data_test.loc[data_test['Screening Status'] == '0', 'Screening Status'] = 'negative'
    data_test.loc[data_test['Screening Status'] == '1', 'Screening Status'] = 'CD8'
    data_test.loc[data_test['Screening Status'] == 'unscreened', 'Screening Status'] = 'not_tested'

    gartner_data = pd.concat([data_test, data_train], ignore_index=True, join='inner')
    gartner_data = gartner_data.astype({'ID': str, 'Mut Epitope': str})

    gartner_patients = set(gartner_data['ID'].unique())

    result_df = pd.DataFrame()
    mutation_data = DataManager.load_filter_data(peptide_type='mutation', dataset='NCI', ml_row_selection=False)

    inhouse_info = []
    for patient in mutation_data['patient'].unique():
        if patient not in gartner_patients:
            continue

        data_p = DataManager.load_filter_data(peptide_type='mutation', patient=patient, ml_row_selection=False)

        gartner_data_sel = gartner_data.loc[gartner_data['ID'] == patient, ]

        gartner_data_sel_SNV = gartner_data_sel.loc[gartner_data_sel['Mutation type'] == 'nonsynonymous SNV', :]
        data_p_SNV = data_p.loc[data_p['mutation_type'] == 'SNV', :]

        d = {'Patient': patient,
             'NCI_train_mut_seq_all in-house': data_p.shape[0],
             'NCI_train_mut_seq_tested in-house': sum(np.logical_or(data_p['response_type'] == 'negative',
                                                      data_p.apply(lambda r: 'CD8' in r['response_type'], axis=1))),
             'NCI_train_mut_seq_imm in-house SNV':
                 sum(np.logical_and(data_p['mutation_type'] == 'SNV',
                                    data_p.apply(lambda r: 'CD8' in r['response_type'], axis=1))),
             'NCI_train_mut_seq_imm in-house': sum(data_p.apply(lambda r: 'CD8' in r['response_type'], axis=1)),
             'NCI_train_mut_seq in-house SNV': sum(data_p['mutation_type'] == 'SNV'),
             'NCI_train_mut_seq in-house InDel': sum(np.logical_or(data_p['mutation_type'] == 'INSERTION',
                                                     data_p['mutation_type'] == 'DELETION')),
             'NCI_train_mut_seq in-house FS': sum(data_p['mutation_type'] == 'FSS'),
             'NCI_train_mut_seq_all Gartner': gartner_data_sel.shape[0],
             'NCI_train_mut_seq_tested Gartner': sum(np.logical_or(data_p['response_type'] == 'negative',
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
             'NCI_train_mut_seq_all in-house/Gartner overlap': get_overlapping_mutation_count(gartner_data_sel, data_p),
             'NCI_train_mut_seq in-house/Gartner SNV overlap': get_overlapping_mutation_count(gartner_data_sel_SNV, data_p_SNV)}

        result_df = pd.concat([result_df, pd.Series(d).to_frame().T], ignore_index=True, axis=0)

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

    result_df.to_csv(os.path.join(GlobalParameters.data_dir, 'Compare_in-house_Gartner_Counts.txt'),
                     header=True, index=False, sep="\t")
