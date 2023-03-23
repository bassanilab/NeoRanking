import argparse

from DataWrangling.DataLoader import *
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Print dataset statistics')

parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='Response types included for testing')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='Mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='Immunogenic response_types included')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')

args = parser.parse_args()

response_types = ['CD8', 'negative', 'not_tested', 'total']


def count_alleles(alleles_list):
    alleles = set()
    for alist in alleles_list:
        alleles = alleles.union(set(alist.split(',')))

    return len(alleles)


def update_allele_counts(allele_dict, row):
    for a in str(row['mutant_best_alleles']).split(','):
        rt = row['response_type']
        allele_dict['Neo-pep'][a] = allele_dict['Neo-pep'][a]+1
        if rt =='CD8':
            allele_dict['Neo-pep_imm'][a] = allele_dict['Neo-pep_imm'][a]+1
            allele_dict['Neo-pep_tested'][a] = allele_dict['Neo-pep_tested'][a]+1
        elif rt == 'negative':
            allele_dict['Neo-pep_non-imm'][a] = allele_dict['Neo-pep_non-imm'][a]+1
            allele_dict['Neo-pep_tested'][a] = allele_dict['Neo-pep_tested'][a]+1
        else:
            allele_dict['Neo-pep_not-tested'][a] = allele_dict['Neo-pep_not-tested'][a]+1


def format_alleles(alleles):
    fas = set()
    for a in alleles:
        fas.add(re.sub("[\\*\\:]", "", a))
    return fas


def count_allele_peptides(data_frame, data_manager):
    alleles = set()
    for p in patients:
        alleles = alleles.union(data_manager.get_classI_allotypes(p))
    alleles = format_alleles(alleles)
    d = {'Neo-pep': dict.fromkeys(alleles, 0), 'Neo-pep_imm': dict.fromkeys(alleles, 0),
         'Neo-pep_non-imm': dict.fromkeys(alleles, 0), 'Neo-pep_not-tested': dict.fromkeys(alleles, 0),
         'Neo-pep_tested': dict.fromkeys(alleles, 0)}
    data_frame.apply(lambda row: update_allele_counts(d, row), axis=1)

    return d


def update_cancer_counts(cancer_dict, row):
    rt = row['response_type']
    cancer_dict[row['Cancer_Type']]['Mut-seq'] = cancer_dict[row['Cancer_Type']]['Mut-seq']+1
    if rt =='CD8':
        cancer_dict[row['Cancer_Type']]['Mut-seq_imm'] = cancer_dict[row['Cancer_Type']]['Mut-seq_imm']+1
        cancer_dict[row['Cancer_Type']]['Mut-seq_tested'] = cancer_dict[row['Cancer_Type']]['Mut-seq_tested']+1
    elif rt == 'negative':
        cancer_dict[row['Cancer_Type']]['Mut-seq_non-imm'] = cancer_dict[row['Cancer_Type']]['Mut-seq_non-imm']+1
        cancer_dict[row['Cancer_Type']]['Mut-seq_tested'] = cancer_dict[row['Cancer_Type']]['Mut-seq_tested']+1
    else:
        cancer_dict[row['Cancer_Type']]['Mut-seq_not-tested'] = cancer_dict[row['Cancer_Type']]['Mut-seq_not-tested']+1


def update_patients(cancer_dict, row):
    cancer_dict[row['Cancer_Type']].add(row['patient'])


def list_cancer_patients(data_frame):
    cancers = data_frame['Cancer_Type'].unique()
    d = dict.fromkeys(cancers, None)
    for c in d:
        d[c] = set()
    data_frame.apply(lambda row: update_patients(d, row), axis=1)

    return d


def count_cancer_peptides(data_frame):
    cancers = data_frame['Cancer_Type'].unique()
    d = dict.fromkeys(cancers, None)
    for c in d:
        d[c] = {'Mut-seq': 0, 'Mut-seq_imm': 0, 'Mut-seq_non-imm': 0, 'Mut-seq_not-tested': 0, 'Mut-seq_tested': 0}
    data_frame.apply(lambda row: update_cancer_counts(d, row), axis=1)

    return d


def get_cancer_statistics(dataset, data_frame):
    list_patients = list_cancer_patients(data_frame)
    cnts_peptides = count_cancer_peptides(data_frame)
    cancer_types = []
    CD8_counts = []
    neg_counts = []
    nt_counts = []
    t_counts = []
    tot_counts = []
    patient_counts = []
    patient_set = []
    for c in cnts_peptides:
        cancer_types.append(c)
        CD8_counts.append(cnts_peptides[c]['Mut-seq_imm'])
        neg_counts.append(cnts_peptides[c]['Mut-seq_non-imm'])
        nt_counts.append(cnts_peptides[c]['Mut-seq_not-tested'])
        t_counts.append(cnts_peptides[c]['Mut-seq_tested'])
        tot_counts.append(cnts_peptides[c]['Mut-seq'])
        patient_counts.append(len(list_patients[c]))
        patient_set.append(list_patients[c])

    d = {'Dataset': dataset, 'Cancer type': cancer_types, 'Patient count': patient_counts, 'Mut-seq count': tot_counts,
         'Mut-seq_imm count': CD8_counts, 'Mut-seq_non-imm count': neg_counts, 'Mut-seq_tested count': t_counts,
         'Mut-seq_not-tested count': nt_counts, 'Patients': patient_set}

    return pd.DataFrame(d)


def count_CD8_peptides_per_mutation(data_frame):
    df = data_frame.loc[:, ['mut_seqid', 'response']]
    return df.groupby('mut_seqid').aggregate('sum')['response'].mean(skipna=True)


def count_CD8_peptides_per_CD8_mutation(data_frame):
    df = data_frame.loc[data_frame['response'] == 1, ['mut_seqid', 'response']]
    return df.groupby('mut_seqid').aggregate('sum')['response'].mean(skipna=True)


def count_CD8_peptides_per_screened_mutation(data_frame):
    idx = data_frame.apply(lambda row: row['response_type'] == 'negative' or row['response_type'] == 'CD8', axis=1)
    df = data_frame.loc[idx, ['mut_seqid', 'response']]
    return df.groupby('mut_seqid').aggregate('sum')['response'].mean(skipna=True)


def get_patient_statistics(patient, patient_group, ml_group, data_frame, peptide_type):
    df = data_frame[data_frame['patient'] == patient]
    cancer_type = data_frame.loc[data_frame['patient'] == patient, 'Cancer_Type'].unique()

    if peptide_type == 'long':
        d = {'Patient': patient,
             'Patient group': patient_group,
             'Cancer type': cancer_type,
             'ML group': ml_group,
             'Mut-seq count': df.shape[0],
             'Mut-seq_imm count': sum(df['response'] == 1),
             'Mut-seq_non-imm count': sum(df['response_type'] == 'negative'),
             'Mut-seq_tested count': sum(df['response_type'] != 'not_tested'),
             'Mut-seq_not-tested count': sum(df['response_type'] == 'not_tested'),
             'Mean mut-seq RNAseq coverage': df['rnaseq_alt_support'].mean(skipna=True),
             'Mean mut-seq_imm RNAseq coverage': df.loc[df['response'] == 1, 'rnaseq_alt_support'].mean(skipna=True),
             'Mean mut-seq_non-imm RNAseq coverage': df.loc[df['response_type'] == 'negative', 'rnaseq_alt_support'].mean(skipna=True),
             'Mean mut-seq_tested RNAseq coverage': df.loc[df['response_type'] != 'not_tested', 'rnaseq_alt_support'].mean(skipna=True),
             'Mean mut-seq_not-tested RNAseq coverage': df.loc[df['response_type'] == 'not_tested', 'rnaseq_alt_support'].mean(skipna=True),
             'Mean mut-seq RNAseq TPM': df['rnaseq_TPM'].mean(skipna=True),
             'Mean mut-seq_imm RNAseq TPM': df.loc[df['response'] == 1, 'rnaseq_TPM'].mean(skipna=True),
             'Mean mut-seq_non-imm RNAseq TPM': df.loc[df['response_type'] == 'negative', 'rnaseq_TPM'].mean(skipna=True),
             'Mean mut-seq_tested RNAseq TPM': df.loc[df['response_type'] != 'not_tested', 'rnaseq_TPM'].mean(skipna=True),
             'Mean mut-seq_not-tested RNAseq TPM': df.loc[df['response_type'] == 'not_tested', 'rnaseq_TPM'].mean(skipna=True),
             }

    else:
        d = {'Patient': patient,
             'Patient group': patient_group,
             'Cancer type': cancer_type,
             'ML group': ml_group,
             'Neo-pep count': df.shape[0],
             'Neo-pep_imm count': sum(df['response'] == 1),
             'Neo-pep_non-imm count': sum(df['response_type'] == 'negative'),
             'Neo-pep_tested count': sum(df['response_type'] != 'not_tested'),
             'Neo-pep_not-tested count': sum(df['response_type'] == 'not_tested'),
             'Mean mut-seq RNAseq coverage': df['rnaseq_alt_support'].mean(skipna=True),
             'Mean mut-seq_imm RNAseq coverage': df.loc[df['response'] == 1, 'rnaseq_alt_support'].mean(skipna=True),
             'Mean mut-seq_non-imm RNAseq coverage': df.loc[df['response_type'] == 'negative', 'rnaseq_alt_support'].mean(skipna=True),
             'Mean mut-seq_tested RNAseq coverage': df.loc[df['response_type'] != 'not_tested', 'rnaseq_alt_support'].mean(skipna=True),
             'Mean mut-seq_not-tested RNAseq coverage': df.loc[df['response_type'] == 'not_tested', 'rnaseq_alt_support'].mean(skipna=True),
             'Mean mut-seq RNAseq TPM': df['rnaseq_TPM'].mean(skipna=True),
             'Mean mut-seq_imm RNAseq TPM': df.loc[df['response'] == 1, 'rnaseq_TPM'].mean(skipna=True),
             'Mean mut-seq_non-imm RNAseq TPM': df.loc[df['response_type'] == 'negative', 'rnaseq_TPM'].mean(skipna=True),
             'Mean mut-seq_tested RNAseq TPM': df.loc[df['response_type'] != 'not_tested', 'rnaseq_TPM'].mean(skipna=True),
             'Mean mut-seq_not-tested RNAseq TPM': df.loc[df['response_type'] == 'not_tested', 'rnaseq_TPM'].mean(skipna=True),
             'Neo-pep allele count': count_alleles(df['mutant_best_alleles']),
             'Neo-pep_imm allele count': count_alleles(df.loc[df['response'] == 1, 'mutant_best_alleles']),
             'Neo-pep_non-imm allele count': count_alleles(df.loc[df['response_type'] == 'negative', 'mutant_best_alleles']),
             'Neo-pep_tested allele count': count_alleles(df.loc[df['response_type'] != 'not_tested', 'mutant_best_alleles']),
             'Neo-pep_not-tested allele count': count_alleles(df.loc[df['response_type'] == 'not_tested', 'mutant_best_alleles']),
             'Mean neo-pep MixMHC rank': df['mutant_rank'].mean(skipna=True),
             'Mean neo-pep_imm MixMHC rank': df.loc[df['response'] == 1, 'mutant_rank'].mean(skipna=True),
             'Mean neo-pep_non-imm MixMHC rank': df.loc[df['response_type'] == 'negative', 'mutant_rank'].mean(skipna=True),
             'Mean neo-pep_tested MixMHC rank': df.loc[df['response_type'] != 'not_tested', 'mutant_rank'].mean(skipna=True),
             'Mean neo-pep_not-tested MixMHC rank': df.loc[df['response_type'] == 'not_tested', 'mutant_rank'].mean(skipna=True),
             'Mean neo-pep_imm count per mut-seq': count_CD8_peptides_per_mutation(df),
             'Mean neo-pep_imm count per mut-seq_imm': count_CD8_peptides_per_CD8_mutation(df),
             'Mean neo-pep_imm count per mut-seq_tested': count_CD8_peptides_per_screened_mutation(df)
             }

    return pd.Series(d)


data_loader = DataLoader(mutation_types=args.mutation_types, response_types=args.response_types,
                         immunogenic=args.immunogenic, min_nr_immuno=0)

pa_output_file = os.path.join(Parameters().get_plot_dir(), "Patient_statistics_{0}.txt".format(args.peptide_type))
open(pa_output_file, 'w').close()
ds_output_file = os.path.join(Parameters().get_plot_dir(), "Dataset_statistics_{0}.txt".format(args.peptide_type))
open(ds_output_file, 'w').close()
ca_output_file = os.path.join(Parameters().get_plot_dir(), "Cancer_statistics_{0}.txt".format(args.peptide_type))
open(ca_output_file, 'w').close()
cnt_output_file = os.path.join(Parameters().get_plot_dir(), "Dataset_counts_{0}.txt".format(args.peptide_type))
open(cnt_output_file, 'w').close()

if args.peptide_type == 'short':
    al_output_file = os.path.join(Parameters().get_plot_dir(), "Allele_statistics_{0}.txt".format(args.peptide_type))
    with open(al_output_file, 'w') as file:
        file.write("Dataset\tAllele\tResponse type\tNeo-pep count\n")


pa_stats_df = pd.DataFrame()
ds_stats_df = pd.DataFrame()
ca_stats_df = pd.DataFrame()

mgr = DataManager()

#for patient_group in ['HiTIDE']:
for patient_group in ['HiTIDE', 'TESLA', 'Rosenberg', 'Gartner_test', 'Gartner_train']:

    patients = get_valid_patients(dataset=patient_group, peptide_type=args.peptide_type)

    data, X, y = data_loader.load_patients(patients, args.input_file_tag, args.peptide_type, verbose=True)

    ds_pa_stats_df = pd.DataFrame()
    for patient in data['patient'].unique():
        ml_group = get_ml_group(patient, args.peptide_type)
        ds_pa_stats_df = ds_pa_stats_df.append(get_patient_statistics(patient, patient_group, ml_group, data,
                                                                      args.peptide_type), ignore_index=True)

    stats = ds_pa_stats_df.describe().reset_index()
    stats.insert(loc=0, column='Dataset', value=patient_group)
    stats = stats.rename(columns={'index': 'Statistics'})
    stats = stats.loc[stats['Statistics'] != 'count']
    ds_stats_df = ds_stats_df.append(stats, ignore_index=True)

    if patient_group not in ['Gartner_test', 'Gartner_train']:
        pa_stats_df = pa_stats_df.append(ds_pa_stats_df, ignore_index=True)

    ca_stats_df = ca_stats_df.append(get_cancer_statistics(patient_group, data), ignore_index=True)

    if args.peptide_type == 'short':
        with open(al_output_file, 'a') as file:
            al_cnts_peptides = count_allele_peptides(data, mgr)
            for rt in al_cnts_peptides:
                for a in al_cnts_peptides[rt]:
                    file.write("{0}\t{1}\t{2}\t{3}\n".format(patient_group, a, rt, al_cnts_peptides[rt][a]))

if args.peptide_type == 'long':
    col_order = ['Patient', 'Cancer type', 'Patient group', 'ML group',
                 'Mut-seq count', 'Mut-seq_imm count', 'Mut-seq_non-imm count', 'Mut-seq_tested count',
                 'Mut-seq_not-tested count', 'Mean mut-seq RNAseq TPM', 'Mean mut-seq_imm RNAseq TPM',
                 'Mean mut-seq_non-imm RNAseq TPM', 'Mean mut-seq_tested RNAseq TPM',
                 'Mean mut-seq_not-tested RNAseq TPM', 'Mean mut-seq RNAseq coverage',
                 'Mean mut-seq_imm RNAseq coverage',  'Mean mut-seq_non-imm RNAseq coverage',
                 'Mean mut-seq_tested RNAseq coverage', 'Mean mut-seq_not-tested RNAseq coverage']

    counts_df = pa_stats_df[['Patient', 'Cancer type', 'Patient group', 'ML group', 'Mut-seq count', 'Mut-seq_imm count',
                             'Mut-seq_non-imm count', 'Mut-seq_tested count', 'Mut-seq_not-tested count']]

else:
    col_order = ['Patient', 'Cancer type', 'Patient group', 'ML group',
                 'Neo-pep count', 'Neo-pep_imm count', 'Neo-pep_non-imm count', 'Neo-pep_tested count',
                 'Neo-pep_not-tested count', 'Mean neo-pep_imm count per mut-seq',
                 'Mean neo-pep_imm count per mut-seq_imm', 'Mean neo-pep_imm count per mut-seq_tested',
                 'Mean mut-seq RNAseq TPM', 'Mean mut-seq_imm RNAseq TPM', 'Mean mut-seq_non-imm RNAseq TPM',
                 'Mean mut-seq_tested RNAseq TPM', 'Mean mut-seq_not-tested RNAseq TPM', 'Mean mut-seq RNAseq coverage',
                 'Mean mut-seq_imm RNAseq coverage', 'Mean mut-seq_non-imm RNAseq coverage',
                 'Mean mut-seq_tested RNAseq coverage', 'Mean mut-seq_not-tested RNAseq coverage',
                 'Neo-pep allele count', 'Neo-pep_imm allele count', 'Neo-pep_non-imm allele count',
                 'Neo-pep_tested allele count', 'Neo-pep_not-tested allele count', 'Mean neo-pep MixMHC rank',
                 'Mean neo-pep_imm MixMHC rank', 'Mean neo-pep_non-imm MixMHC rank',
                 'Mean neo-pep_tested MixMHC rank', 'Mean neo-pep_not-tested MixMHC rank']

    counts_df = pa_stats_df[['Patient', 'Cancer type', 'Patient group', 'ML group', 'Neo-pep count', 'Neo-pep_imm count',
                             'Neo-pep_non-imm count', 'Neo-pep_tested count', 'Neo-pep_not-tested count']]

pa_stats_df = pa_stats_df[col_order]
pa_stats_df.to_csv(pa_output_file, sep="\t", index=False, header=True)
ds_stats_df.to_csv(ds_output_file, sep="\t", index=False, header=True)
ca_stats_df.to_csv(ca_output_file, sep="\t", index=False, header=True)

counts_df = counts_df.groupby(['Patient group', 'ML group']).sum()
counts_df.to_csv(cnt_output_file, sep="\t", index=True, header=True)



