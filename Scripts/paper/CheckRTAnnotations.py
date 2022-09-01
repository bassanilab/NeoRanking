import argparse
import pandas as pd
import time

from DataWrangling.DataLoader import *
from Utils.Util_fct import *


def match_long_short_rows(row, mutant_seq):
    mut_pos = row['pep_mut_start']
    mut_seq_long = row['mutant_seq']
    start = max(mut_pos-len(mutant_seq), 0)
    end = min(mut_pos+len(mutant_seq), len(mut_seq_long))
    return mutant_seq in mut_seq_long[start:end]


def get_long_rt(row, df_long, response_types_dict):
    mut_seqid = row['mut_seqid']
    if mut_seqid in response_types_dict:
        df_long_sel = response_types_dict[mut_seqid]
    else:
        df_long_sel = \
            df_long.loc[(df_long['mutant_id'] == mut_seqid), :]
        response_types_dict[mut_seqid] = df_long_sel
        # if df_long_sel.shape[0] > 1:
        #     print(df_long_sel.shape[0])

    idx = df_long_sel.apply(match_long_short_rows, args=(row['mutant_seq'], ), axis=1)
    response_types = df_long_sel.loc[idx, 'response_type']

    if any(response_types == 'CD8'):
        return 'CD8'
    elif any(response_types == 'negative'):
        return 'negative'
    else:
        return 'not_tested'


def convert_peptide_id_long(row):
    fields = row['peptide_id'].split('|')
    return fields[0]+":"+fields[2]


def add_long_rt_to_short(data_long, data_short):
    assert 'peptide_id' in data_long.columns, "No peptide_id in data"
    
    data_long.loc[:, 'mutant_id'] = data_long.apply(convert_peptide_id_long, axis=1)

    data_long = data_long[['peptide_id', 'mutant_id', 'pep_mut_start', 'mutant_seq', 'response_type']]
    data_short = data_short[['peptide_id', 'mut_seqid', 'mutant_seq', 'response_type']]

    response_types = {}
    long_rt = data_short.apply(get_long_rt, args=(data_long, response_types), axis=1)

    data_short.loc[:, 'long_response_type'] = long_rt

    return data_short


def main(args):

    patients_long = \
        get_valid_patients(patients=args.patients, peptide_type='long') \
            if args.patients and len(args.patients) > 0 else get_valid_patients(peptide_type='long')

    patients_short = \
        get_valid_patients(patients=args.patients, peptide_type='short') \
            if args.patients and len(args.patients) > 0 else get_valid_patients(peptide_type='short')

    nr_included = len(patients_short.intersection(patients_long))

    data_loader = DataLoader(mutation_types=args.mutation_types, response_types=args.response_types,
                             immunogenic=args.immunogenic, min_nr_immuno=0)

    ct_all = pd.DataFrame()
    data_all = None
    for i, p in enumerate(patients_short):
        data_long, X_long, y_long = \
            data_loader.load_patients(p, args.input_file_tag, 'long', verbose=True)
        data_short, X_short, y_short = \
            data_loader.load_patients(p, args.input_file_tag, 'short', verbose=True)

        start_time = time.time()
        data_short = add_long_rt_to_short(data_long, data_short)
        print("--- %s seconds ---" % (time.time() - start_time))

        print("Processed patient {0} ({1} of {2}) in {3}secs".format(p, i+1, len(patients_short),
                                                                     time.time() - start_time))
        ct = pd.crosstab(index=data_short['response_type'], columns=data_short['long_response_type'])
        if data_all is None:
            ct_all = ct
            data_all = data_short
        else:
            data_all = data_all.append(data_short, ignore_index=True)
            ct_all = ct.add(ct_all, fill_value=0)

        print(ct)

    print("nr long: {0}, nr short: {1}, overlap: {2}".format(len(patients_long), len(patients_short), nr_included))

    print(ct_all)

    data_all.to_csv(path_or_buf=args.data_file, sep="\t", index=False, header=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Write all patient data into one file')
    parser.add_argument('-f', '--data_file', type=str, help='Tabular output file')
    parser.add_argument('-p', '--patients', type=str, nargs='+', help='Patients')
    parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included for testing')
    parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
    parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
    parser.add_argument('-i', '--input_file_tag', type=str, default='rt',
                        help='File tag for neodisc input file (patient)_(peptide_type)_(input_file_tag).txt')

    args = parser.parse_args()

    for arg in vars(args):
        print(arg, getattr(args, arg))

    main(args)
