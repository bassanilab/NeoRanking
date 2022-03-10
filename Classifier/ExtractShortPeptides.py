from Classify_functions import *
from Utils.Parameters import *

cl_vars = ['Peptide_Group', 'Filter', 'Samples', 'Nb_Samples', 'peptide_id','chromosome', 'genomic_coord', 'ref', 'alt',
           'gene', 'aa_wt', 'protein_coord', 'aa_mutant', 'Sample_Tissue_expression_GTEx', 'rnaseq_TPM', 'mutant_seq',
           'mutant_best_alleles', 'mutant_rank', 'wt_seq', 'wt_rank', 'Final_Peptide_Class', 'database_entry',
           'GENOMIC_MUTATION_ID (COSMIC)', 'FATHMM prediction (COSMIC)', 'gene_driver_Intogen',
           'mutation_driver_statement_Intogen']


def add_best_n_peptides(data, i, nr_peptides, short_peptides_df, j, cl_vars, max_netmhc_rank):

    idx = data.index[i]
    nr_brest_peptides = nr_peptides[i]

    best_short_peptides = get_best_peptides(data, idx, nr_brest_peptides, max_netmhc_rank)

    for k in np.arange(best_short_peptides.shape[0]):
        for cl_var in cl_vars:
            if cl_var in data.columns:
                short_peptides_df.loc[j+k, cl_var] = data.loc[idx, cl_var]

        for c in best_short_peptides.columns:
            short_peptides_df.loc[j+k, c] = best_short_peptides.loc[best_short_peptides.index[k], c]

        short_peptides_df.loc[j+k, 'Peptide_Group'] = "NA"
        short_peptides_df.loc[j+k, 'Final_Peptide_Class'] = "HLAI"

    return j + best_short_peptides.shape[0]


def get_best_peptides(data, i, nr_brest_peptides, max_netmhc_rank):

    columns = ['position', 'mutant_seq', 'mutant_best_alleles', 'mutant_rank', 'wt_seq', 'wt_best_alleles', 'wt_rank']
    short_peptide_res = pd.DataFrame(columns=columns, index=np.arange(nr_brest_peptides))
    k = 0
    for j in np.arange(max_netmhc_rank):
        p = str(data.loc[i, 'mut_peptide_'+str(j)])
        if 9 <= len(p) <= 11:
            short_peptide_res.loc[k, 'position'] = data.loc[i, 'mut_peptide_pos_'+str(j)]
            short_peptide_res.loc[k, 'mutant_seq'] = p
            short_peptide_res.loc[k, 'mutant_best_alleles'] = data.loc[i, 'mut_allele_'+str(j)]
            short_peptide_res.loc[k, 'mutant_rank'] = \
                0.5*(data.loc[i, 'mut_Rank_EL_'+str(j)]+data.loc[i, 'mut_Rank_BA_'+str(j)])
            short_peptide_res.loc[k, 'wt_seq'] = data.loc[i, 'wt_peptide_'+str(j)]
            short_peptide_res.loc[k, 'wt_best_alleles'] = data.loc[i, 'wt_allele_'+str(j)]
            short_peptide_res.loc[k, 'wt_rank'] = \
                0.5*(data.loc[i, 'wt_Rank_EL_'+str(j)]+data.loc[i, 'wt_Rank_BA_'+str(j)])
            k += 1

    short_peptide_res.dropna(axis=0, inplace=True)
    short_peptide_res.sort_values(by=['mutant_rank'], ascending=True, inplace=True)

    keep_index = []
    used_positions = []
    used_alleles = []
    for j in short_peptide_res.index:
        pos = int(short_peptide_res.loc[j, 'position'])
        allele = str(short_peptide_res.loc[j, 'mutant_best_alleles'])
        new_pos = all([abs(p - pos) > 2 for p in used_positions])
        new_allele = all([allele != a for a in used_alleles])
        ratio = short_peptide_res.loc[j, 'mutant_rank']/short_peptide_res.loc[j, 'wt_rank']
        same_allele = allele == short_peptide_res.loc[j, 'wt_best_alleles']
        mut_wt_diff = (not same_allele and ratio < 1.5) or ratio < 2.5
        mut_rank = short_peptide_res.loc[j, 'mutant_rank'] < 3

        if (new_pos or new_allele) and mut_wt_diff and mut_rank:
            keep_index.append(j)
            used_alleles.append(allele)
            used_positions.append(pos)

    short_peptide_res = short_peptide_res.loc[keep_index, :]
    short_peptide_res.sort_values(by=['mutant_rank'], ascending=True, inplace=True)

    columns = ['mutant_seq', 'mutant_best_alleles', 'mutant_rank', 'wt_seq', 'wt_rank']
    short_peptide_res = short_peptide_res[columns]

    return short_peptide_res


def calc_nr_peptides(scores, tot_nr_peptides, max_nr_peptides):

    peptide_cnts = np.divide(scores, max(scores)/(max_nr_peptides-1))
    peptide_cnts = np.array(np.round(peptide_cnts), dtype='int32')
    peptide_cnts = np.add(peptide_cnts, 1)

    tot_cnt = 0
    for i in np.arange(len(peptide_cnts)):
        avail_peptide_cnt = tot_nr_peptides - tot_cnt
        if avail_peptide_cnt < 0:
            if peptide_cnts[i-1] > 0:
               peptide_cnts[i-1] += avail_peptide_cnt
            peptide_cnts[i] = 0
        tot_cnt += peptide_cnts[i]

    return peptide_cnts


def update_peptide_cnts(i, spare_cnts, peptide_cnts, max_peptides_per_mut):
    peptide_cnts[i] -= spare_cnts
    j = i+1
    while spare_cnts > 0 and j < len(peptide_cnts):
        if peptide_cnts[j] < max_peptides_per_mut:
            peptide_cnts[j] += 1
            spare_cnts -= 1
        j += 1


def extract_short_peptides_all(data, max_nr):

    max_peptides_per_mut = 4
    netMHCpan_ranks = [int(c[c.rfind('_')+1:]) for c in data.columns if 'mut_peptide_pos_' in c]
    peptide_cnts = calc_nr_peptides(data['ML_pred'], max_nr, min(len(netMHCpan_ranks), max_peptides_per_mut))

    short_peptide_res = pd.DataFrame(columns=cl_vars, index=np.arange(max_nr))

    j = 0
    i = 0
    while j < max_nr and i < data.shape[0]:
        j_prev = j
        j = add_best_n_peptides(data, i, peptide_cnts, short_peptide_res, j, cl_vars, len(netMHCpan_ranks))

        spare_cnts = peptide_cnts[i] - j + j_prev
        if spare_cnts > 0:
            update_peptide_cnts(i, spare_cnts, peptide_cnts, max_peptides_per_mut)
        i += 1

    return short_peptide_res


file_dir = '/home/localadmin/Priorization/results/'
file_tag = 'netmhc_stab_chop_ML_SVM_13'

patients = ['13LN', '1HU3', '1G4T', '1I3M', '0YM1n']
for patient in patients:

    if Parameters().patient_exist(patient):
        print("extract short peptides for: "+str(patient))
        df, y = load_data(patient, file_tag)

        if df is None:
            continue

        short_peptide_df = extract_short_peptides_all(df, 120)

        out_file = os.path.join(Parameters().get_result_dir(), patient+"_"+file_tag+"_short_peptides.txt")
        short_peptide_df.to_csv(path_or_buf=out_file, sep="\t", index=False, header=True)

    else:
        print("Patient "+str(patient)+" does not have a data file. Skip.")

