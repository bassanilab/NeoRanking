import argparse
import pandas as pd
import re

from Utils.Util_fct import *
from Utils.DataManager import DataManager
from DataWrangling.DataTransformer import DataTransformer


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


def count_allele_peptides(neopep_data_):
    alleles = set()
    for p in DataManager.allotypes.Patient:
        alleles = alleles.union(DataManager.get_classI_allotypes(p))
    alleles = format_alleles(alleles)
    allele_cnts = {'Neo-pep': dict.fromkeys(alleles, 0),
                   'Neo-pep_imm': dict.fromkeys(alleles, 0),
                   'Neo-pep_non-imm': dict.fromkeys(alleles, 0),
                   'Neo-pep_not-tested': dict.fromkeys(alleles, 0),
                   'Neo-pep_tested': dict.fromkeys(alleles, 0)}
    neopep_data_.apply(lambda row: update_allele_counts(allele_cnts, row), axis=1)

    return allele_cnts


def count_CD8_peptides_per_mutation(neopep_data_):
    neopep_grouped = neopep_data_.groupby(['chromosome', 'genomic_coord', 'ref', 'alt'])
    cnts = []
    for key, grp in neopep_grouped:
        cnts.append(grp['response_type'].apply(lambda e: int(e == 'CD8')).sum())

    if len(cnts) > 0:
        return np.mean(cnts)
    else:
        return np.nan


def count_CD8_peptides_per_CD8_mutation(neopep_data_, mutation_data_):
    mutation_data_ = mutation_data_.loc[mutation_data_['response_type'] == 'CD8', :]
    mutation_grouped = mutation_data_.groupby(['chromosome', 'genomic_coord', 'ref', 'alt'])
    neopep_grouped = neopep_data_.groupby(['chromosome', 'genomic_coord', 'ref', 'alt'])
    cnts = []
    for key, g in mutation_grouped:
        if key not in neopep_grouped.groups.keys():
            continue
        grp = neopep_grouped.get_group(key)
        cnts.append(grp['response_type'].apply(lambda e: int(e == 'CD8')).sum())

    if len(cnts) > 0:
        return np.mean(cnts)
    else:
        return np.nan


def count_CD8_peptides_per_screened_mutation(neopep_data_, mutation_data_):
    mutation_data_ = mutation_data_.loc[mutation_data_['response_type'] != 'not_tested', :]
    mutation_grouped = mutation_data_.groupby(['chromosome', 'genomic_coord', 'ref', 'alt'])
    neopep_grouped = neopep_data_.groupby(['chromosome', 'genomic_coord', 'ref', 'alt'])
    cnts = []
    for key, g in mutation_grouped:
        if key not in neopep_grouped.groups.keys():
            continue
        grp = neopep_grouped.get_group(key)
        cnts.append(grp['response_type'].apply(lambda e: int(e == 'CD8')).sum())

    if len(cnts) > 0:
        return np.mean(cnts)
    else:
        return np.nan


def get_patient_mutation_statistics(patient_, dataset_, mutation_data_):
    mutation_data_ = mutation_data_[mutation_data_['patient'] == patient_]
    cancer_type = mutation_data_.loc[mutation_data_['patient'] == patient_, 'Cancer_Type'].unique()
    ml_group = mutation_data_.loc[mutation_data_['patient'] == patient_, 'train_test'].iloc[0]

    mutation_stats = {
        'Patient': patient_,
        'Patient group': dataset_,
        'Cancer type': cancer_type,
        'ML group': ml_group,
        'Mut-seq count': mutation_data_.shape[0],
        'Mut-seq_imm count': sum(mutation_data_['response_type'] == 'CD8'),
        'Mut-seq_non-imm count': sum(mutation_data_['response_type'] == 'negative'),
        'Mut-seq_tested count': sum(mutation_data_['response_type'] != 'not_tested'),
        'Mut-seq_not-tested count': sum(mutation_data_['response_type'] == 'not_tested'),
        'Mean mut-seq RNAseq coverage': mutation_data_['rnaseq_alt_support'].mean(skipna=True),
        'Mean mut-seq_imm RNAseq coverage':
            mutation_data_.loc[mutation_data_['response_type'] == 'CD8', 'rnaseq_alt_support'].mean(skipna=True),
        'Mean mut-seq_non-imm RNAseq coverage':
            mutation_data_.loc[mutation_data_['response_type'] == 'negative', 'rnaseq_alt_support'].mean(skipna=True),
        'Mean mut-seq_tested RNAseq coverage':
            mutation_data_.loc[mutation_data_['response_type'] != 'not_tested', 'rnaseq_alt_support'].mean(skipna=True),
        'Mean mut-seq_not-tested RNAseq coverage':
            mutation_data_.loc[mutation_data_['response_type'] == 'not_tested', 'rnaseq_alt_support'].mean(skipna=True),
        'Mean mut-seq RNAseq TPM': mutation_data_['rnaseq_TPM'].mean(skipna=True),
        'Mean mut-seq_imm RNAseq TPM':
            mutation_data_.loc[mutation_data_['response_type'] == 'CD8', 'rnaseq_TPM'].mean(skipna=True),
        'Mean mut-seq_non-imm RNAseq TPM':
            mutation_data_.loc[mutation_data_['response_type'] == 'negative', 'rnaseq_TPM'].mean(skipna=True),
        'Mean mut-seq_tested RNAseq TPM':
            mutation_data_.loc[mutation_data_['response_type'] != 'not_tested', 'rnaseq_TPM'].mean(skipna=True),
        'Mean mut-seq_not-tested RNAseq TPM':
            mutation_data_.loc[mutation_data_['response_type'] == 'not_tested', 'rnaseq_TPM'].mean(skipna=True),
    }

    return pd.Series(mutation_stats).to_frame().T


def get_patient_neopep_statistics(patient_, dataset_, neopep_data_, mutation_data_):
    neopep_data_ = neopep_data_[neopep_data_['patient'] == patient_]
    mutation_data_ = mutation_data_[mutation_data_['patient'] == patient_]
    cancer_type = mutation_data_.loc[mutation_data_['patient'] == patient_, 'Cancer_Type'].unique()
    ml_group = mutation_data_.loc[mutation_data_['patient'] == patient_, 'train_test'].iloc[0]

    neopep_stats = {
        'Patient': patient_,
        'Patient group': dataset_,
        'Cancer type': cancer_type,
        'ML group': ml_group,
        'Neo-pep count': neopep_data_.shape[0],
        'Neo-pep_imm count': sum(neopep_data_['response_type'] == 'CD8'),
        'Neo-pep_non-imm count': sum(neopep_data_['response_type'] == 'negative'),
        'Neo-pep_tested count': sum(neopep_data_['response_type'] != 'not_tested'),
        'Neo-pep_not-tested count': sum(neopep_data_['response_type'] == 'not_tested'),
        'Mean mut-seq RNAseq coverage': neopep_data_['rnaseq_alt_support'].mean(skipna=True),
        'Mean mut-seq_imm RNAseq coverage':
            neopep_data_.loc[neopep_data_['response_type'] == 'CD8', 'rnaseq_alt_support'].mean(skipna=True),
        'Mean mut-seq_non-imm RNAseq coverage':
            neopep_data_.loc[neopep_data_['response_type'] == 'negative', 'rnaseq_alt_support'].mean(skipna=True),
        'Mean mut-seq_tested RNAseq coverage':
            neopep_data_.loc[neopep_data_['response_type'] != 'not_tested', 'rnaseq_alt_support'].mean(skipna=True),
        'Mean mut-seq_not-tested RNAseq coverage':
            neopep_data_.loc[neopep_data_['response_type'] == 'not_tested', 'rnaseq_alt_support'].mean(skipna=True),
        'Mean mut-seq RNAseq TPM': neopep_data_['rnaseq_TPM'].mean(skipna=True),
        'Mean mut-seq_imm RNAseq TPM':
            neopep_data_.loc[neopep_data_['response_type'] == 'CD8', 'rnaseq_TPM'].mean(skipna=True),
        'Mean mut-seq_non-imm RNAseq TPM':
            neopep_data_.loc[neopep_data_['response_type'] == 'negative', 'rnaseq_TPM'].mean(skipna=True),
        'Mean mut-seq_tested RNAseq TPM':
            neopep_data_.loc[neopep_data_['response_type'] != 'not_tested', 'rnaseq_TPM'].mean(skipna=True),
        'Mean mut-seq_not-tested RNAseq TPM':
            neopep_data_.loc[neopep_data_['response_type'] == 'not_tested', 'rnaseq_TPM'].mean(skipna=True),
        'Neo-pep allele count': count_alleles(neopep_data_['mutant_best_alleles']),
        'Neo-pep_imm allele count':
            count_alleles(neopep_data_.loc[neopep_data_['response_type'] == 'CD8', 'mutant_best_alleles']),
        'Neo-pep_non-imm allele count':
            count_alleles(neopep_data_.loc[neopep_data_['response_type'] == 'negative', 'mutant_best_alleles']),
        'Neo-pep_tested allele count':
            count_alleles(neopep_data_.loc[neopep_data_['response_type'] != 'not_tested', 'mutant_best_alleles']),
        'Neo-pep_not-tested allele count':
            count_alleles(neopep_data_.loc[neopep_data_['response_type'] == 'not_tested', 'mutant_best_alleles']),
        'Mean neo-pep MixMHC rank': neopep_data_['mutant_rank'].mean(skipna=True),
        'Mean neo-pep_imm MixMHC rank':
            neopep_data_.loc[neopep_data_['response_type'] == 'CD8', 'mutant_rank'].mean(skipna=True),
        'Mean neo-pep_non-imm MixMHC rank':
            neopep_data_.loc[neopep_data_['response_type'] == 'negative', 'mutant_rank'].mean(skipna=True),
        'Mean neo-pep_tested MixMHC rank':
            neopep_data_.loc[neopep_data_['response_type'] != 'not_tested', 'mutant_rank'].mean(skipna=True),
        'Mean neo-pep_not-tested MixMHC rank':
            neopep_data_.loc[neopep_data_['response_type'] == 'not_tested', 'mutant_rank'].mean(skipna=True),
        'Mean neo-pep_imm count per mut-seq': count_CD8_peptides_per_mutation(neopep_data_),
        'Mean neo-pep_imm count per mut-seq_imm': count_CD8_peptides_per_CD8_mutation(neopep_data_, mutation_data_),
        'Mean neo-pep_imm count per mut-seq_tested': count_CD8_peptides_per_screened_mutation(neopep_data_, mutation_data_)
    }

    return pd.Series(neopep_stats).to_frame().T


if __name__ == "__main__":

    pa_neopep_output_file = os.path.join(GlobalParameters.plot_dir, "Patient_statistics_neopep.txt")
    open(pa_neopep_output_file, 'w').close()
    cnt_neopep_output_file = os.path.join(GlobalParameters.plot_dir, "Dataset_counts_neopep.txt")
    open(cnt_neopep_output_file, 'w').close()
    pa_mutation_output_file = os.path.join(GlobalParameters.plot_dir, "Patient_statistics_mutation.txt")
    open(pa_mutation_output_file, 'w').close()
    cnt_mutation_output_file = os.path.join(GlobalParameters.plot_dir, "Dataset_counts_mutation.txt")
    open(cnt_mutation_output_file, 'w').close()

    al_output_file = os.path.join(GlobalParameters.plot_dir, "Allele_statistics_neopep.txt")
    with open(al_output_file, 'w') as file:
        file.write("Dataset\tAllele\tResponse type\tNeo-pep count\n")

    pa_mutation_stats_df = pd.DataFrame()
    pa_neopep_stats_df = pd.DataFrame()

    for dataset in ['HiTIDE', 'TESLA', 'NCI']:

        neopep_data = DataManager.filter_original_data(peptide_type='neopep', dataset=dataset, ml_row_selection=True)
        mutation_data = DataManager.filter_original_data(peptide_type='mutation', dataset=dataset, ml_row_selection=True)

        ds_pa_mutation_stats_df = pd.DataFrame()
        for patient in mutation_data['patient'].unique():
            ds_pa_mutation_stats_df = \
                pd.concat([ds_pa_mutation_stats_df, get_patient_mutation_statistics(patient, dataset, mutation_data)],
                          ignore_index=True, axis=0)

        ds_pa_neopep_stats_df = pd.DataFrame()
        for patient in neopep_data['patient'].unique():
            ds_pa_neopep_stats_df = \
                pd.concat([ds_pa_neopep_stats_df, get_patient_neopep_statistics(patient, dataset, neopep_data, mutation_data)],
                          ignore_index=True, axis=0)

        pa_mutation_stats_df = pd.concat([pa_mutation_stats_df, ds_pa_mutation_stats_df], ignore_index=True)
        pa_neopep_stats_df = pd.concat([pa_neopep_stats_df, ds_pa_neopep_stats_df], ignore_index=True)

        with open(al_output_file, 'a') as file:
            al_cnts_peptides = count_allele_peptides(neopep_data)
            for rt in al_cnts_peptides:
                for a in al_cnts_peptides[rt]:
                    file.write("{0}\t{1}\t{2}\t{3}\n".format(dataset, a, rt, al_cnts_peptides[rt][a]))

    col_order = ['Patient', 'Cancer type', 'Patient group', 'ML group',
                 'Mut-seq count', 'Mut-seq_imm count', 'Mut-seq_non-imm count', 'Mut-seq_tested count',
                 'Mut-seq_not-tested count', 'Mean mut-seq RNAseq TPM', 'Mean mut-seq_imm RNAseq TPM',
                 'Mean mut-seq_non-imm RNAseq TPM', 'Mean mut-seq_tested RNAseq TPM',
                 'Mean mut-seq_not-tested RNAseq TPM', 'Mean mut-seq RNAseq coverage',
                 'Mean mut-seq_imm RNAseq coverage',  'Mean mut-seq_non-imm RNAseq coverage',
                 'Mean mut-seq_tested RNAseq coverage', 'Mean mut-seq_not-tested RNAseq coverage']
    pa_mutation_stats_df = pa_mutation_stats_df[col_order]

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
    pa_neopep_stats_df = pa_neopep_stats_df[col_order]

    pa_mutation_stats_df.to_csv(pa_mutation_output_file, sep="\t", index=False, header=True)
    pa_neopep_stats_df.to_csv(pa_neopep_output_file, sep="\t", index=False, header=True)



