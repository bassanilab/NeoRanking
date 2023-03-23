import os.path


class Parameters:

    def __init__(self):
        self.base_dir = "/home/localadmin/Priorization"
        self.data_dir = os.path.join(self.base_dir, "data")
        self.exe_dir = "/home/localadmin/Programmes"
        self.result_dir = os.path.join(self.base_dir, "results")
#        self.result_dir = os.path.join(self.base_dir, "results_ipMSDB_test")
        self.plot_dir = os.path.join(self.base_dir, "Plots")
        self.pickle_dir = os.path.join(self.base_dir, "Classifiers")

        self.data_immunogenicity_info_file = \
            os.path.join(self.data_dir, "Data_immunogenicity_info.txt")
        self.data_validity_file = os.path.join(self.data_dir, "Data_validity_info.txt")
        self.allo_file = os.path.join(self.data_dir, "hla", "HLA_allotypes.txt")
        self.cat_to_num_info_files = {'short': {}, 'long': {}}
        self.cat_to_num_info_files['short']['NCI_train'] = os.path.join(self.data_dir, "Cat_to_num_info_short_NCI_train.txt")
        self.cat_to_num_info_files['long']['NCI_train'] = os.path.join(self.data_dir, "Cat_to_num_info_long_NCI_train.txt")
        self.cat_to_num_info_files['short']['NCI'] = os.path.join(self.data_dir, "Cat_to_num_info_short_NCI_all.txt")
        self.cat_to_num_info_files['long']['NCI'] = os.path.join(self.data_dir, "Cat_to_num_info_long_NCI_all.txt")
        self.protein_seq_file_37 = os.path.join(self.data_dir, "fasta", "gencode.v38lift37.pc_translations.reformatted.fa")
        self.protein_seq_file_38 = os.path.join(self.data_dir, "fasta", "Homo_sapiens.GRCh38.pep.all.fa")
        self.htide_info_file = \
            os.path.join(self.data_dir, "immunogenicity", "AgDisc_peptides_2023-02-21_all.txt")
        # self.htide_info_file = \
        #     os.path.join(self.data_dir, "immunogenicity", "20220110_Immunogenicity_testing_database.xlsx")
        self.gartner_info_long_train = os.path.join(self.data_dir, "immunogenicity", "NmersTrainingSet.txt")
        self.gartner_info_long_test = os.path.join(self.data_dir, "immunogenicity", "NmersTestingSet.txt")
        self.gartner_info_short_train = os.path.join(self.data_dir, "immunogenicity", "MmpsTrainingSet.txt")
        self.gartner_info_short_test = os.path.join(self.data_dir, "immunogenicity", "MmpsTestingSet.txt")
        self.gartner_info_file = os.path.join(self.data_dir, "immunogenicity", "43018_2021_197_MOESM2_ESM.xlsx")
        self.parkhurst_info_file_long = \
            os.path.join(self.data_dir, "immunogenicity", "214162_2_supp_5551286_ps8cdj.xlsx")
        self.tesla_info_short_1 = os.path.join(self.data_dir, "immunogenicity", "mmc4.xlsx")
        self.tesla_info_short_2 = os.path.join(self.data_dir, "immunogenicity", "mmc7.xlsx")
        self.tesla_results = os.path.join(self.data_dir, "misc", "mmc5.xlsx")
        self.gartner_long_results = os.path.join(self.data_dir, "immunogenicity", "Gartner_nmers_ranking.txt")

        self.ipmsdb_file = os.path.join(self.data_dir, "misc", "Comet_peptides_standard_protein_summary.txt")
        self.human_fasta_file = os.path.join(self.data_dir, "misc", "uniprot-human.fasta")
        self.virus_fasta_file = os.path.join(self.data_dir, "misc", "virus.proteome.fasta")
        self.human_allele_score_dist_file = os.path.join(self.data_dir, "hla", "allele_score_data_human.txt")
        self.virus_allele_score_dist_file = os.path.join(self.data_dir, "hla", "allele_score_data_virus.txt")

        self.aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

        self.env = \
            {'NETMHCpan': '/home/localadmin/Programmes/netMHCpan-4.1/Linux_x86_64',
             'NETMHCstabpan': '/home/localadmin/Programmes/netMHCstabpan-1.0/Linux_x86_64',
             'NETCHOP': '/home/localadmin/Programmes/netchop-3.1/Linux_x86_64',
             'NETCTLpan': '/home/localadmin/Programmes/netCTLpan-1.1/Linux_x86_64',
             'TMPDIR': '/home/localadmin/tmp'}

        self.features = \
            ['peptide_id', 'mut_seqid', 'TumorContent',  'callers neoDisc_identification', 'ref', 'alt', 'mutation_type',
             'Clonality', 'Zygosity', 'VAF', 'gene', 'database_entry', 'mutation', 'Sample_Tissue_expression_GTEx',
             'Cancer_Type', 'TCGA_Cancer_expression', 'aa_wt', 'protein_coord', 'aa_mutant', 'rnaseq_TPM', 'CCF',
             'rnaseq_gene_expression_quartile', 'rnaseq_coverage', 'rnaseq_detected_variants', 'rnaseq_ref_support',
             'rnaseq_alt_support', 'pep_mut_start', 'pep_mut_end', 'seq_len', 'mutant_seq', 'wt_seq', 'Nb_Samples',
             'MIN_MUT_RANK_CI_MIXMHC', 'COUNT_MUT_RANK_CI_MIXMHC', 'WT_RANK_CI_MIXMHC', 'WT_BEST_RANK_CI_MIXMHC',
             'MIN_MUT_RANK_CI_PRIME', 'COUNT_MUT_RANK_CI_PRIME', 'WT_RANK_CI_PRIME', 'WT_BEST_RANK_CI_PRIME',
             'MIN_MUT_RANK_CI_netMHCpan', 'COUNT_MUT_RANK_CI_netMHCpan', 'WT_RANK_CI_netMHCpan',
             'WT_BEST_RANK_CI_netMHCpan',  'MIN_MUT_RANK_CII', 'INTRACELLULAR_LOCATIONS', 'EXTRACELLULAR_LOCATIONS',
             'CSCAPE_score', 'Rank_mutation_frequency (COSMIC)', 'GENOMIC_MUTATION_ID (COSMIC)',
             'LEGACY_MUTATION_ID (COSMIC)', 'FATHMM prediction (COSMIC)', 'FATHMM score (COSMIC)',
             'Frequency_(counts/nb_samples) (COSMIC)', 'gene_driver_Intogen', 'nb_mutations_in_gene_Intogen',
             'nb_mutations_same_position_Intogen', 'nb_same_mutation_Intogen', 'mutation_driver_statement_Intogen',
             'GTEx_all_tissues_expression_median', 'GTEx_all_tissues_expression_mean', 'bestWTMatchScore_I',
             'bestWTMatchOverlap_I', 'bestMutationScore_I', 'bestWTPeptideCount_I', 'bestWTMatchType_I',
             'response_type', 'MIN_MUT_RANK_CII', 'WT_BEST_RANK_CII', 'COUNT_MUT_RANK_CII',
             'mut_peptide_pos_0', 'mut_allele_0', 'mut_peptide_0', 'mut_core_0',
             'mut_Of_0', 'mut_Gp_0', 'mut_Gl_0', 'mut_Ip_0', 'mut_Il_0', 'mut_Icore_0',
             'mut_Score_EL_0', 'mut_Rank_EL_0', 'mut_Score_BA_0', 'mut_Rank_BA_0', 'mut_ic50_0',
             'mut_pos_in_peptide_0', 'mut_nr_strong_binders_0', 'mut_nr_weak_binders_0',
             'mut_nr_weak_binding_alleles_0', 'wt_allele_0', 'wt_peptide_0', 'wt_Score_EL_0', 'wt_Rank_EL_0',
             'wt_Score_BA_0', 'wt_Rank_BA_0', 'wt_ic50_0', 'mut_peptide_pos_1', 'mut_allele_1', 'mut_peptide_1',
             'mut_core_1', 'mut_Of_1', 'mut_Gp_1', 'mut_Gl_1', 'mut_Ip_1', 'mut_Il_1', 'mut_Icore_1',
             'mut_Score_EL_1', 'mut_Rank_EL_1', 'mut_Score_BA_1', 'mut_Rank_BA_1', 'mut_ic50_1',
             'mut_pos_in_peptide_1', 'wt_allele_1', 'wt_peptide_1', 'wt_Score_EL_1', 'wt_Rank_EL_1',
             'wt_Score_BA_1', 'wt_Rank_BA_1', 'wt_ic50_1', 'mut_peptide_pos_2', 'mut_allele_2', 'mut_peptide_2',
             'mut_core_2', 'mut_Of_2', 'mut_Gp_2', 'mut_Gl_2', 'mut_Ip_2', 'mut_Il_2', 'mut_Icore_2',
             'mut_Score_EL_2', 'mut_Rank_EL_2', 'mut_Score_BA_2', 'mut_Rank_BA_2', 'mut_ic50_2',
             'mut_pos_in_peptide_2', 'wt_allele_2', 'wt_peptide_2', 'wt_Score_EL_2', 'wt_Rank_EL_2',
             'wt_Score_BA_2', 'wt_Rank_BA_2', 'wt_ic50_2', 'mut_peptide_pos_3', 'mut_allele_3', 'mut_peptide_3',
             'mut_core_3', 'mut_Of_3', 'mut_Gp_3', 'mut_Gl_3', 'mut_Ip_3', 'mut_Il_3', 'mut_Icore_3',
             'mut_Score_EL_3', 'mut_Rank_EL_3', 'mut_Score_BA_3', 'mut_Rank_BA_3', 'mut_ic50_3',
             'mut_pos_in_peptide_3', 'wt_allele_3', 'wt_peptide_3', 'wt_Score_EL_3', 'wt_Rank_EL_3',
             'wt_Score_BA_3', 'wt_Rank_BA_3', 'wt_ic50_3', 'mut_peptide_pos_4', 'mut_allele_4', 'mut_peptide_4',
             'mut_core_4', 'mut_Of_4', 'mut_Gp_4', 'mut_Gl_4', 'mut_Ip_4', 'mut_Il_4', 'mut_Icore_4',
             'mut_Score_EL_4', 'mut_Rank_EL_4', 'mut_Score_BA_4', 'mut_Rank_BA_4', 'mut_ic50_4',
             'mut_pos_in_peptide_4', 'wt_allele_4', 'wt_peptide_4', 'wt_Score_EL_4', 'wt_Rank_EL_4',
             'wt_Score_BA_4', 'wt_Rank_BA_4', 'wt_ic50_4', 'mut_Stab_Score_0', 'mut_Thalf_0', 'mut_Rank_Stab_0',
             'wt_Stab_Score_0', 'wt_Thalf_0', 'wt_Rank_Stab_0', 'mut_Stab_Score_1', 'mut_Thalf_1', 'mut_Rank_Stab_1',
             'wt_Stab_Score_1', 'wt_Thalf_1', 'wt_Rank_Stab_1', 'mut_Stab_Score_2', 'mut_Thalf_2', 'mut_Rank_Stab_2',
             'wt_Stab_Score_2', 'wt_Thalf_2', 'wt_Rank_Stab_2', 'mut_Stab_Score_3', 'mut_Thalf_3', 'mut_Rank_Stab_3',
             'wt_Stab_Score_3', 'wt_Thalf_3', 'wt_Rank_Stab_3', 'mut_Stab_Score_4', 'mut_Thalf_4', 'mut_Rank_Stab_4',
             'wt_Stab_Score_4', 'wt_Thalf_4', 'wt_Rank_Stab_4', 'wt_netchop_Ct_score_0', 'mut_netchop_Ct_score_0',
             'wt_netchop_Nt_score_0', 'mut_netchop_Nt_score_0', 'wt_netchop_Int_score_0', 'mut_netchop_Int_score_0',
             'wt_netchop_Ct_score_1', 'mut_netchop_Ct_score_1', 'wt_netchop_Nt_score_1', 'mut_netchop_Nt_score_1',
             'wt_netchop_Int_score_1', 'mut_netchop_Int_score_1', 'wt_netchop_Ct_score_2', 'mut_netchop_Ct_score_2',
             'wt_netchop_Nt_score_2', 'mut_netchop_Nt_score_2', 'wt_netchop_Int_score_2', 'mut_netchop_Int_score_2',
             'wt_netchop_Ct_score_3', 'mut_netchop_Ct_score_3', 'wt_netchop_Nt_score_3', 'mut_netchop_Nt_score_3',
             'wt_netchop_Int_score_3', 'mut_netchop_Int_score_3', 'wt_netchop_Ct_score_4', 'mut_netchop_Ct_score_4',
             'wt_netchop_Nt_score_4', 'mut_netchop_Nt_score_4', 'wt_netchop_Int_score_4', 'mut_netchop_Int_score_4',
             'mut_TAP_score_0', 'wt_TAP_score_0', 'mut_TAP_score_1', 'wt_TAP_score_1', 'mut_TAP_score_2',
             'wt_TAP_score_2', 'mut_TAP_score_3', 'wt_TAP_score_3', 'mut_TAP_score_4', 'wt_TAP_score_4',
             'mut_is_binding_pos_0', 'wt_binding_score_0', 'mut_binding_score_0', 'mut_is_binding_pos_1',
             'wt_binding_score_1', 'mut_binding_score_1', 'mut_is_binding_pos_2', 'wt_binding_score_2',
             'mut_binding_score_2', 'mut_is_binding_pos_3', 'wt_binding_score_3', 'mut_binding_score_3',
             'mut_is_binding_pos_4', 'wt_binding_score_4', 'mut_binding_score_4',
             'mut_aa_coeff_0', 'wt_aa_coeff_0', 'next_best_EL_mut_ranks', 'next_best_BA_mut_ranks',
             'mut_aa_coeff_1', 'wt_aa_coeff_1', 'mut_aa_coeff_2', 'wt_aa_coeff_2',
             'mut_aa_coeff_3', 'wt_aa_coeff_3', 'mut_aa_coeff_4', 'wt_aa_coeff_4',
             'DAI_0', 'DAI_1', 'DAI_2', 'DAI_3', 'DAI_4', 'mut_allele_propensity_0', 'mut_allele_propensity_1',
             'mut_allele_propensity_2', 'mut_allele_propensity_3', 'mut_allele_propensity_4', 'wt_allele_propensity_0',
             'wt_allele_propensity_1', 'wt_allele_propensity_2', 'wt_allele_propensity_3', 'wt_allele_propensity_4',
             'mutant_rank', 'wt_best_rank', 'mutant_best_alleles', 'mutant_other_significant_alleles',
             'number_overlaping_HLA_II', 'mutant_best_alleles_PRIME', 'mutant_rank_PRIME', 'wt_best_rank_PRIME',
             'mutant_best_alleles_netMHCpan', 'mutant_rank_netMHCpan', 'wt_best_rank_netMHCpan', 'Sample_Tissue', 'DAI',
             'mut_Stab_Score', 'mut_Thalf', 'mut_Rank_Stab', 'wt_Stab_Score', 'wt_Thalf', 'wt_Rank_Stab',
             'mut_netchop_score', 'mut_netchop_score_ct', 'mut_netchop_score_nt', 'mut_netchop_score_int',
             'mut_is_binding_pos', 'mut_binding_score', 'TAP_score', 'mut_aa_coeff', 'wt_aa_coeff',
             'mut_allele_propensity', 'rank_in_mutation', 'number_included_HLA_I', 'peptide_score', 'mutation_score',
             'mutation_rank', 'DAI_NetMHC', 'DAI_MixMHC', 'DAI_NetStab', 'DAI_MixMHC_mbp', 'INCLUDED_SHORT_PEPTIDES',
             'DAI_aa_coeff']

        self.ml_features = \
            ['TumorContent', 'CCF', 'Clonality', 'Nb_Samples',
             'Zygosity', 'VAF', 'mutation', 'Sample_Tissue_expression_GTEx', 'Cancer_Type',
             'TCGA_Cancer_expression', 'aa_wt', 'protein_coord', 'aa_mutant', 'rnaseq_TPM',
             'rnaseq_gene_expression_quartile', 'rnaseq_coverage', 'rnaseq_detected_variants', 'rnaseq_ref_support',
             'rnaseq_alt_support', 'pep_mut_start', 'pep_mut_end', 'seq_len',
             'MIN_MUT_RANK_CI_MIXMHC', 'COUNT_MUT_RANK_CI_MIXMHC', 'WT_RANK_CI_MIXMHC', 'WT_BEST_RANK_CI_MIXMHC',
             'MIN_MUT_RANK_CI_PRIME', 'COUNT_MUT_RANK_CI_PRIME', 'WT_RANK_CI_PRIME', 'WT_BEST_RANK_CI_PRIME',
             'MIN_MUT_RANK_CI_netMHCpan', 'COUNT_MUT_RANK_CI_netMHCpan', 'WT_RANK_CI_netMHCpan',
             'WT_BEST_RANK_CI_netMHCpan',  'MIN_MUT_RANK_CII', 'INTRACELLULAR_LOCATIONS',
             'CSCAPE_score', 'Rank_mutation_frequency (COSMIC)', 'FATHMM prediction (COSMIC)', 'FATHMM score (COSMIC)',
             'Frequency_(counts/nb_samples) (COSMIC)', 'gene_driver_Intogen', 'nb_mutations_in_gene_Intogen',
             'nb_mutations_same_position_Intogen', 'nb_same_mutation_Intogen', 'mutation_driver_statement_Intogen',
             'GTEx_all_tissues_expression_median', 'GTEx_all_tissues_expression_mean', 'bestWTMatchScore_I',
             'bestWTMatchOverlap_I', 'bestMutationScore_I', 'bestWTPeptideCount_I', 'bestWTMatchType_I',
             'response_type', 'MIN_MUT_RANK_CII', 'WT_BEST_RANK_CII', 'COUNT_MUT_RANK_CII',
             'mut_peptide_pos_0', 'mut_allele_0', 'mut_peptide_0', 'mut_core_0',
             'mut_Of_0', 'mut_Gp_0', 'mut_Gl_0', 'mut_Ip_0', 'mut_Il_0', 'mut_Icore_0',
             'mut_Score_EL_0', 'mut_Rank_EL_0', 'mut_Score_BA_0', 'mut_Rank_BA_0', 'mut_ic50_0',
             'mut_pos_in_peptide_0', 'mut_nr_strong_binders_0', 'mut_nr_weak_binders_0',
             'mut_nr_weak_binding_alleles_0', 'wt_allele_0', 'wt_peptide_0', 'wt_Score_EL_0', 'wt_Rank_EL_0',
             'wt_Score_BA_0', 'wt_Rank_BA_0', 'wt_ic50_0', 'mut_peptide_pos_1', 'mut_allele_1', 'mut_peptide_1',
             'mut_core_1', 'mut_Of_1', 'mut_Gp_1', 'mut_Gl_1', 'mut_Ip_1', 'mut_Il_1', 'mut_Icore_1',
             'mut_Score_EL_1', 'mut_Rank_EL_1', 'mut_Score_BA_1', 'mut_Rank_BA_1', 'mut_ic50_1',
             'mut_pos_in_peptide_1', 'wt_allele_1', 'wt_peptide_1', 'wt_Score_EL_1', 'wt_Rank_EL_1',
             'wt_Score_BA_1', 'wt_Rank_BA_1', 'wt_ic50_1', 'mut_peptide_pos_2', 'mut_allele_2', 'mut_peptide_2',
             'mut_core_2', 'mut_Of_2', 'mut_Gp_2', 'mut_Gl_2', 'mut_Ip_2', 'mut_Il_2', 'mut_Icore_2',
             'mut_Score_EL_2', 'mut_Rank_EL_2', 'mut_Score_BA_2', 'mut_Rank_BA_2', 'mut_ic50_2',
             'mut_pos_in_peptide_2', 'wt_allele_2', 'wt_peptide_2', 'wt_Score_EL_2', 'wt_Rank_EL_2',
             'wt_Score_BA_2', 'wt_Rank_BA_2', 'wt_ic50_2', 'mut_peptide_pos_3', 'mut_allele_3', 'mut_peptide_3',
             'mut_core_3', 'mut_Of_3', 'mut_Gp_3', 'mut_Gl_3', 'mut_Ip_3', 'mut_Il_3', 'mut_Icore_3',
             'mut_Score_EL_3', 'mut_Rank_EL_3', 'mut_Score_BA_3', 'mut_Rank_BA_3', 'mut_ic50_3',
             'mut_pos_in_peptide_3', 'wt_allele_3', 'wt_peptide_3', 'wt_Score_EL_3', 'wt_Rank_EL_3',
             'wt_Score_BA_3', 'wt_Rank_BA_3', 'wt_ic50_3', 'mut_peptide_pos_4', 'mut_allele_4', 'mut_peptide_4',
             'mut_core_4', 'mut_Of_4', 'mut_Gp_4', 'mut_Gl_4', 'mut_Ip_4', 'mut_Il_4', 'mut_Icore_4',
             'mut_Score_EL_4', 'mut_Rank_EL_4', 'mut_Score_BA_4', 'mut_Rank_BA_4', 'mut_ic50_4',
             'mut_pos_in_peptide_4', 'wt_allele_4', 'wt_peptide_4', 'wt_Score_EL_4', 'wt_Rank_EL_4',
             'wt_Score_BA_4', 'wt_Rank_BA_4', 'wt_ic50_4', 'mut_Stab_Score_0', 'mut_Thalf_0', 'mut_Rank_Stab_0',
             'wt_Stab_Score_0', 'wt_Thalf_0', 'wt_Rank_Stab_0', 'mut_Stab_Score_1', 'mut_Thalf_1', 'mut_Rank_Stab_1',
             'wt_Stab_Score_1', 'wt_Thalf_1', 'wt_Rank_Stab_1', 'mut_Stab_Score_2', 'mut_Thalf_2', 'mut_Rank_Stab_2',
             'wt_Stab_Score_2', 'wt_Thalf_2', 'wt_Rank_Stab_2', 'mut_Stab_Score_3', 'mut_Thalf_3', 'mut_Rank_Stab_3',
             'wt_Stab_Score_3', 'wt_Thalf_3', 'wt_Rank_Stab_3', 'mut_Stab_Score_4', 'mut_Thalf_4', 'mut_Rank_Stab_4',
             'wt_Stab_Score_4', 'wt_Thalf_4', 'wt_Rank_Stab_4', 'wt_netchop_Ct_score_0', 'mut_netchop_Ct_score_0',
             'wt_netchop_Nt_score_0', 'mut_netchop_Nt_score_0', 'wt_netchop_Int_score_0', 'mut_netchop_Int_score_0',
             'wt_netchop_Ct_score_1', 'mut_netchop_Ct_score_1', 'wt_netchop_Nt_score_1', 'mut_netchop_Nt_score_1',
             'wt_netchop_Int_score_1', 'mut_netchop_Int_score_1', 'wt_netchop_Ct_score_2', 'mut_netchop_Ct_score_2',
             'wt_netchop_Nt_score_2', 'mut_netchop_Nt_score_2', 'wt_netchop_Int_score_2', 'mut_netchop_Int_score_2',
             'wt_netchop_Ct_score_3', 'mut_netchop_Ct_score_3', 'wt_netchop_Nt_score_3', 'mut_netchop_Nt_score_3',
             'wt_netchop_Int_score_3', 'mut_netchop_Int_score_3', 'wt_netchop_Ct_score_4', 'mut_netchop_Ct_score_4',
             'wt_netchop_Nt_score_4', 'mut_netchop_Nt_score_4', 'wt_netchop_Int_score_4', 'mut_netchop_Int_score_4',
             'mut_TAP_score_0', 'wt_TAP_score_0', 'mut_TAP_score_1', 'wt_TAP_score_1', 'mut_TAP_score_2',
             'wt_TAP_score_2', 'mut_TAP_score_3', 'wt_TAP_score_3', 'mut_TAP_score_4', 'wt_TAP_score_4',
             'mut_is_binding_pos_0', 'wt_binding_score_0', 'mut_binding_score_0', 'mut_is_binding_pos_1',
             'wt_binding_score_1', 'mut_binding_score_1', 'mut_is_binding_pos_2', 'wt_binding_score_2',
             'mut_binding_score_2', 'mut_is_binding_pos_3', 'wt_binding_score_3', 'mut_binding_score_3',
             'mut_is_binding_pos_4', 'wt_binding_score_4', 'mut_binding_score_4',
             'mut_aa_coeff_0', 'wt_aa_coeff_0', 'next_best_EL_mut_ranks', 'next_best_BA_mut_ranks',
             'mut_aa_coeff_1', 'wt_aa_coeff_1', 'mut_aa_coeff_2', 'wt_aa_coeff_2',
             'mut_aa_coeff_3', 'wt_aa_coeff_3', 'mut_aa_coeff_4', 'wt_aa_coeff_4',
             'DAI_0', 'DAI_1', 'DAI_2', 'DAI_3', 'DAI_4', 'mut_allele_propensity_0', 'mut_allele_propensity_1',
             'mut_allele_propensity_2', 'mut_allele_propensity_3', 'mut_allele_propensity_4', 'wt_allele_propensity_0',
             'wt_allele_propensity_1', 'wt_allele_propensity_2', 'wt_allele_propensity_3', 'wt_allele_propensity_4',
             'mutant_rank', 'wt_best_rank', 'mutant_best_alleles', 'mutant_other_significant_alleles',
             'number_overlaping_HLA_II', 'mutant_best_alleles_PRIME', 'mutant_rank_PRIME', 'wt_best_rank_PRIME',
             'mutant_best_alleles_netMHCpan', 'mutant_rank_netMHCpan', 'wt_best_rank_netMHCpan', 'Sample_Tissue', 'DAI',
             'mut_Stab_Score', 'mut_Thalf', 'mut_Rank_Stab', 'wt_Stab_Score', 'wt_Thalf', 'wt_Rank_Stab',
             'mut_netchop_score', 'mut_netchop_score_ct', 'mut_netchop_score_nt', 'mut_netchop_score_int',
             'mut_is_binding_pos', 'mut_binding_score', 'TAP_score', 'mut_aa_coeff', 'wt_aa_coeff',
             'mut_allele_propensity', 'rank_in_mutation', 'number_included_HLA_I', 'peptide_score', 'mutation_score',
             'mutation_rank', 'DAI_NetMHC', 'DAI_MixMHC', 'DAI_NetStab', 'DAI_MixMHC_mbp', 'bestWTMatchType_I',
             'DAI_aa_coeff']

        self.num_features = \
            ['VAF', 'rnaseq_ref_support', 'rnaseq_alt_support', 'CCF', 'Nb_Samples',
             'Sample_Tissue_expression_GTEx', 'protein_coord', 'TCGA_Cancer_expression', 'rnaseq_TPM',
             'COUNT_MUT_RANK_CI_netMHCpan', 'COUNT_MUT_RANK_CI_MIXMHC', 'WT_RANK_CI_MIXMHC', 'WT_BEST_RANK_CI_MIXMHC',
             'MIN_MUT_RANK_CI_PRIME', 'COUNT_MUT_RANK_CI_PRIME', 'WT_RANK_CI_PRIME', 'WT_BEST_RANK_CI_PRIME',
             'CSCAPE_score', 'Rank_mutation_frequency..COSMIC.', 'FATHMM.score..COSMIC.',
             'nb_mutations_in_gene_Intogen', 'nb_mutations_same_position_Intogen', 'nb_same_mutation_Intogen',
             'Frequency_.counts.nb_samples...COSMIC.', 'GTEx_all_tissues_expression_median',
             'GTEx_all_tissues_expression_mean', 'bestWTMatchScore_I', 'bestWTMatchOverlap_I', 'bestMutationScore_I',
             'bestWTPeptideCount_I', 'MIN_MUT_RANK_CII', 'WT_BEST_RANK_CII', 'COUNT_MUT_RANK_CII',
             'mut_peptide_pos_0', 'mut_Of_0', 'mut_Gp_0', 'mut_Gl_0',
             'mut_Ip_0', 'mut_Il_0', 'mut_Score_EL_0', 'mut_Rank_EL_0', 'mut_Score_BA_0', 'mut_Rank_BA_0', 'mut_ic50_0',
             'mut_pos_in_peptide_0', 'mut_nr_strong_binders_0', 'mut_nr_weak_binders_0',
             'mut_nr_weak_binding_alleles_0', 'wt_Score_EL_0', 'wt_Rank_EL_0', 'wt_Score_BA_0', 'wt_Rank_BA_0',
             'wt_ic50_0', 'mut_peptide_pos_1', 'mut_Of_1', 'mut_Gp_1', 'mut_Gl_1', 'mut_Ip_1', 'mut_Il_1',
             'mut_Score_EL_1', 'mut_Rank_EL_1', 'mut_Score_BA_1', 'mut_Rank_BA_1', 'mut_ic50_1',
             'mut_pos_in_peptide_1', 'wt_Score_EL_1', 'wt_Rank_EL_1', 'wt_Score_BA_1', 'wt_Rank_BA_1',
             'wt_ic50_1', 'mut_peptide_pos_2', 'mut_Of_2', 'mut_Gp_2', 'mut_Gl_2', 'mut_Ip_2', 'mut_Il_2',
             'mut_Score_EL_2', 'mut_Rank_EL_2', 'mut_Score_BA_2', 'mut_Rank_BA_2', 'mut_ic50_2',
             'mut_pos_in_peptide_2', 'wt_Score_EL_2', 'wt_Rank_EL_2', 'wt_Score_BA_2', 'wt_Rank_BA_2',
             'wt_ic50_2', 'mut_peptide_pos_3', 'mut_Of_3', 'mut_Gp_3', 'mut_Gl_3', 'mut_Ip_3', 'mut_Il_3',
             'mut_Score_EL_3', 'mut_Rank_EL_3', 'mut_Score_BA_3', 'mut_Rank_BA_3', 'mut_ic50_3',
             'mut_pos_in_peptide_3', 'wt_Score_EL_3', 'wt_Rank_EL_3', 'wt_Score_BA_3', 'wt_Rank_BA_3',
             'wt_ic50_3', 'mut_peptide_pos_4', 'mut_Of_4', 'mut_Gp_4', 'mut_Gl_4', 'mut_Ip_4', 'mut_Il_4',
             'mut_Score_EL_4', 'mut_Rank_EL_4', 'mut_Score_BA_4', 'mut_Rank_BA_4', 'mut_ic50_4',
             'mut_pos_in_peptide_4', 'wt_Score_EL_4', 'wt_Rank_EL_4', 'wt_Score_BA_4', 'wt_Rank_BA_4',
             'wt_ic50_4', 'mut_Stab_Score_0', 'mut_Thalf_0', 'mut_Rank_Stab_0', 'wt_Stab_Score_0', 'wt_Thalf_0',
             'wt_Rank_Stab_0', 'mut_Stab_Score_1', 'mut_Thalf_1', 'mut_Rank_Stab_1', 'wt_Stab_Score_1', 'wt_Thalf_1',
             'wt_Rank_Stab_1', 'mut_Stab_Score_2', 'mut_Thalf_2', 'mut_Rank_Stab_2', 'wt_Stab_Score_2', 'wt_Thalf_2',
             'wt_Rank_Stab_2', 'mut_Stab_Score_3', 'mut_Thalf_3', 'mut_Rank_Stab_3', 'wt_Stab_Score_3', 'wt_Thalf_3',
             'wt_Rank_Stab_3', 'mut_Stab_Score_4', 'mut_Thalf_4', 'mut_Rank_Stab_4', 'wt_Stab_Score_4', 'wt_Thalf_4',
             'wt_Rank_Stab_4', 'wt_netchop_Ct_score_0', 'mut_netchop_Ct_score_0',
             'wt_netchop_Nt_score_0', 'mut_netchop_Nt_score_0', 'wt_netchop_Int_score_0', 'mut_netchop_Int_score_0',
             'wt_netchop_Ct_score_1', 'mut_netchop_Ct_score_1', 'wt_netchop_Nt_score_1', 'mut_netchop_Nt_score_1',
             'wt_netchop_Int_score_1', 'mut_netchop_Int_score_1', 'wt_netchop_Ct_score_2', 'mut_netchop_Ct_score_2',
             'wt_netchop_Nt_score_2', 'mut_netchop_Nt_score_2', 'wt_netchop_Int_score_2', 'mut_netchop_Int_score_2',
             'wt_netchop_Ct_score_3', 'mut_netchop_Ct_score_3', 'wt_netchop_Nt_score_3', 'mut_netchop_Nt_score_3',
             'wt_netchop_Int_score_3', 'mut_netchop_Int_score_3', 'wt_netchop_Ct_score_4', 'mut_netchop_Ct_score_4',
             'wt_netchop_Nt_score_4', 'mut_netchop_Nt_score_4', 'wt_netchop_Int_score_4', 'mut_netchop_Int_score_4',
             'mut_TAP_score_0', 'wt_TAP_score_0', 'mut_TAP_score_1', 'wt_TAP_score_1', 'mut_TAP_score_2',
             'wt_TAP_score_2', 'mut_TAP_score_3', 'wt_TAP_score_3', 'mut_TAP_score_4', 'wt_TAP_score_4',
             'wt_binding_score_0', 'mut_binding_score_0', 'wt_binding_score_1', 'mut_binding_score_1',
             'wt_binding_score_2', 'mut_binding_score_2', 'wt_binding_score_3', 'mut_binding_score_3',
             'wt_binding_score_4', 'mut_binding_score_4',
             'mut_aa_coeff_0', 'wt_aa_coeff_0', 'next_best_EL_mut_ranks', 'next_best_BA_mut_ranks',
             'DAI_0', 'DAI_1', 'DAI_2', 'DAI_3', 'DAI_4', 'mut_allele_propensity_0', 'mut_allele_propensity_1',
             'mut_allele_propensity_2', 'mut_allele_propensity_3', 'mut_allele_propensity_4', 'wt_allele_propensity_0',
             'wt_allele_propensity_1', 'wt_allele_propensity_2', 'wt_allele_propensity_3', 'wt_allele_propensity_4',
             'mutant_rank', 'wt_best_rank', 'mutant_rank_PRIME', 'wt_best_rank_PRIME', 'mutant_rank_netMHCpan',
             'wt_best_rank_netMHCpan', 'DAI', 'mut_Stab_Score', 'mut_Thalf', 'mut_Rank_Stab', 'wt_Stab_Score',
             'wt_Thalf', 'wt_Rank_Stab', 'mut_netchop_score', 'mut_netchop_score_ct', 'mut_netchop_score_nt',
             'mut_netchop_score_int', 'mut_binding_score', 'TAP_score', 'mut_aa_coeff', 'wt_aa_coeff',
             'mut_allele_propensity', 'peptide_score', 'mutation_score', 'DAI_NetMHC', 'DAI_MixMHC', 'DAI_NetStab',
             'DAI_MixMHC_mbp', 'DAI_aa_coeff']

        # order relation for numerical features (if '>' ('>') missing values are assumed to be low  (high)
        # if relation is '=', a random value is drawn. used only for missing value imputation.
        self.num_features_order = {
            'VAF': '>', 'rnaseq_ref_support': '<', 'rnaseq_alt_support': '>', 'Sample_Tissue_expression_GTEx': '>',
            'TCGA_Cancer_expression': '>', 'rnaseq_TPM': '>', 'CCF': '>', 'seq_len': '=',
            'protein_coord': '=', 'pep_mut_start': '=', 'pep_mut_end': '=', 'MIN_MUT_RANK_CI_MIXMHC': '<',
            'COUNT_MUT_RANK_CI_netMHCpan': '>',
            'COUNT_MUT_RANK_CI_MIXMHC': '>', 'WT_RANK_CI_MIXMHC': '<', 'WT_BEST_RANK_CI_MIXMHC': '<',
            'MIN_MUT_RANK_CI_PRIME': '<', 'COUNT_MUT_RANK_CI_PRIME': '>', 'WT_RANK_CI_PRIME': '<',
            'WT_BEST_RANK_CI_PRIME': '<', 'CSCAPE_score': '>', 'Rank_mutation_frequency..COSMIC.': '>',
            'FATHMM.score..COSMIC.': '>', 'Frequency_.counts.nb_samples...COSMIC.': '>',
            'nb_mutations_in_gene_Intogen': '>', 'nb_mutations_same_position_Intogen': '>',
            'nb_same_mutation_Intogen': '>', 'GTEx_all_tissues_expression_median': '>',
            'GTEx_all_tissues_expression_mean': '>', 'bestWTMatchScore_I': '>', 'bestWTMatchOverlap_I': '>',
            'bestMutationScore_I': '>', 'bestWTPeptideCount_I': '>', 'mut_peptide_pos_0': '=', 'mut_Of_0': '=',
            'mut_Gp_0': '=', 'mut_Gl_0': '=', 'mut_Ip_0': '=', 'mut_Il_0': '=', 'mut_Score_EL_0': '>',
            'mut_Rank_EL_0': '<', 'mut_Score_BA_0': '>', 'mut_Rank_BA_0': '<', 'mut_ic50_0': '<',
            'mut_pos_in_peptide_0': '=', 'mut_nr_strong_binders_0': '>', 'mut_nr_weak_binders_0': '>',
            'mut_nr_weak_binding_alleles_0': '>', 'wt_Score_EL_0': '>', 'wt_Rank_EL_0': '<', 'wt_Score_BA_0': '>',
            'wt_Rank_BA_0': '<', 'wt_ic50_0': '<', 'mut_peptide_pos_1': '=', 'mut_Of_1': '=', 'mut_Gp_1': '=',
            'mut_Gl_1': '=', 'mut_Ip_1': '=', 'mut_Il_1': '=', 'mut_Score_EL_1': '>', 'mut_Rank_EL_1': '<',
            'mut_Score_BA_1': '>', 'mut_Rank_BA_1': '<', 'mut_ic50_1': '<', 'mut_pos_in_peptide_1': '=',
            'mut_nr_strong_binders_1': '>', 'mut_nr_weak_binders_1': '>', 'mut_nr_weak_binding_alleles_1': '>',
            'wt_Score_EL_1': '>', 'wt_Rank_EL_1': '<', 'wt_Score_BA_1': '>', 'wt_Rank_BA_1': '<', 'wt_ic50_1': '<',
            'mut_peptide_pos_2': '=', 'mut_Of_2': '=', 'mut_Gp_2': '=', 'mut_Gl_2': '=', 'mut_Ip_2': '=',
            'mut_Il_2': '=', 'mut_Score_EL_2': '>', 'mut_Rank_EL_2': '<', 'mut_Score_BA_2': '>', 'mut_Rank_BA_2': '<',
            'mut_ic50_2': '<', 'mut_pos_in_peptide_2': '=', 'mut_nr_strong_binders_2': '>',
            'mut_nr_weak_binders_2': '>', 'mut_nr_weak_binding_alleles_2': '>', 'wt_Score_EL_2': '>',
            'wt_Rank_EL_2': '<', 'wt_Score_BA_2': '>', 'wt_Rank_BA_2': '<', 'wt_ic50_2': '<', 'mut_peptide_pos_3': '=',
            'mut_Of_3': '=', 'mut_Gp_3': '=', 'mut_Gl_3': '=', 'mut_Ip_3': '=', 'mut_Il_3': '=',
            'mut_Score_EL_3': '>', 'mut_Rank_EL_3': '<', 'mut_Score_BA_3': '>', 'mut_Rank_BA_3': '<', 'mut_ic50_3': '<',
            'mut_pos_in_peptide_3': '=', 'mut_nr_strong_binders_3': '>', 'mut_nr_weak_binders_3': '>',
            'mut_nr_weak_binding_alleles_3': '>', 'wt_Score_EL_3': '>', 'wt_Rank_EL_3': '<', 'wt_Score_BA_3': '>',
            'wt_Rank_BA_3': '<', 'wt_ic50_3': '<', 'mut_peptide_pos_4': '=', 'mut_Of_4': '=', 'mut_Gp_4': '=',
            'mut_Gl_4': '=', 'mut_Ip_4': '=', 'mut_Il_4': '=', 'mut_Score_EL_4': '>', 'mut_Rank_EL_4': '<',
            'mut_Score_BA_4': '>', 'mut_Rank_BA_4': '<', 'mut_ic50_4': '<', 'mut_pos_in_peptide_4': '=',
            'mut_nr_strong_binders_4': '>', 'mut_nr_weak_binders_4': '>', 'mut_nr_weak_binding_alleles_4': '>',
            'wt_Score_EL_4': '>', 'wt_Rank_EL_4': '<', 'wt_Score_BA_4': '>', 'wt_Rank_BA_4': '<', 'wt_ic50_4': '<',
            'mut_Stab_Score_0': '>', 'mut_Thalf_0': '>', 'mut_Rank_Stab_0': '<', 'wt_Stab_Score_0': '>',
            'wt_Thalf_0': '>',
            'wt_Rank_Stab_0': '<', 'mut_Stab_Score_1': '>', 'mut_Thalf_1': '>', 'mut_Rank_Stab_1': '<',
            'wt_Stab_Score_1': '>', 'wt_Thalf_1': '>', 'wt_Rank_Stab_1': '<', 'mut_Stab_Score_2': '>',
            'mut_Thalf_2': '>', 'mut_Rank_Stab_2': '<', 'wt_Stab_Score_2': '>', 'wt_Thalf_2': '>',
            'wt_Rank_Stab_2': '<',
            'mut_Stab_Score_3': '>', 'mut_Thalf_3': '>', 'mut_Rank_Stab_3': '<', 'wt_Stab_Score_3': '>',
            'wt_Thalf_3': '>',
            'wt_Rank_Stab_3': '<', 'mut_Stab_Score_4': '>', 'mut_Thalf_4': '>', 'mut_Rank_Stab_4': '<',
            'wt_Stab_Score_4': '>', 'wt_Thalf_4': '>', 'wt_Rank_Stab_4': '<', 'wt_netchop_Ct_score_0': '>',
            'mut_netchop_Ct_score_0': '=', 'wt_netchop_Nt_score_0': '=', 'mut_netchop_Nt_score_0': '=',
            'wt_netchop_Int_score_0': '=', 'mut_netchop_Int_score_0': '=', 'wt_netchop_Ct_score_1': '=',
            'mut_netchop_Ct_score_1': '=', 'wt_netchop_Nt_score_1': '=', 'mut_netchop_Nt_score_1': '=',
            'wt_netchop_Int_score_1': '=', 'mut_netchop_Int_score_1': '=', 'wt_netchop_Ct_score_2': '=',
            'mut_netchop_Ct_score_2': '=', 'wt_netchop_Nt_score_2': '=', 'mut_netchop_Nt_score_2': '=',
            'wt_netchop_Int_score_2': '=', 'mut_netchop_Int_score_2': '=', 'wt_netchop_Ct_score_3': '=',
            'mut_netchop_Ct_score_3': '=', 'wt_netchop_Nt_score_3': '=', 'mut_netchop_Nt_score_3': '=',
            'wt_netchop_Int_score_3': '=', 'mut_netchop_Int_score_3': '=', 'wt_netchop_Ct_score_4': '=',
            'mut_netchop_Ct_score_4': '=', 'wt_netchop_Nt_score_4': '=', 'mut_netchop_Nt_score_4': '=',
            'wt_netchop_Int_score_4': '=', 'mut_netchop_Int_score_4': '=', 'mut_TAP_score_0': '>', 'wt_TAP_score_0': '>',
            'mut_TAP_score_1': '>', 'wt_TAP_score_1': '>', 'mut_TAP_score_2': '>', 'wt_TAP_score_2': '>',
            'mut_TAP_score_3': '>', 'wt_TAP_score_3': '>', 'mut_TAP_score_4': '>', 'wt_TAP_score_4': '>',
            'wt_binding_score_0': '=', 'mut_binding_score_0': '=', 'wt_binding_score_1': '=', 'mut_binding_score_1': '=',
            'wt_binding_score_2': '=', 'mut_binding_score_2': '=', 'wt_binding_score_3': '=', 'mut_binding_score_3': '=',
            'wt_binding_score_4': '=', 'mut_binding_score_4': '=', 'mut_aa_coeff_0': '=', 'wt_aa_coeff_0': '=',
            'mut_aa_coeff_1': '=', 'wt_aa_coeff_1': '=', 'mut_aa_coeff_2': '=', 'wt_aa_coeff_2': '=',
            'mut_aa_coeff_3': '=', 'wt_aa_coeff_3': '=',
            'mut_aa_coeff_4': '=', 'wt_aa_coeff_4': '=', 'next_best_EL_mut_ranks': '<', 'next_best_BA_mut_ranks': '<',
            'DAI_0': '=', 'DAI_1': '=', 'DAI_2': '=', 'DAI_3': '=', 'DAI_4': '=', 'mut_allele_propensity_0': '=',
            'mut_allele_propensity_1': '=', 'mut_allele_propensity_2': '=', 'mut_allele_propensity_3': '=',
            'mut_allele_propensity_4': '=', 'wt_allele_propensity_0': '=', 'wt_allele_propensity_1': '=',
            'wt_allele_propensity_2': '=', 'wt_allele_propensity_3': '=', 'wt_allele_propensity_4': '=',
            'rnaseq_gene_expression_quartile': '>', 'mutant_rank': '<', 'wt_best_rank': '<',
            'number_overlaping_HLA_II': '>', 'mutant_rank_PRIME': '<', 'number_included_HLA_I': '>',
            'wt_best_rank_PRIME': '<',  'mutant_rank_netMHCpan': '<',  'wt_best_rank_netMHCpan': '<',
            'DAI': '>', 'mutant_other_significant_alleles': '>', 'mut_Stab_Score': '>', 'mut_Thalf': '>',
            'mut_Rank_Stab': '<', 'wt_Stab_Score': '>', 'wt_Thalf': '>', 'wt_Rank_Stab': '<', 'mut_netchop_score': '>',
            'mut_netchop_score_ct': '>', 'mut_netchop_score_nt': '>', 'mut_netchop_score_int': '<',
            'mut_binding_score': '>', 'TAP_score': '>', 'mut_aa_coeff': '>', 'wt_aa_coeff': '>',
            'mut_allele_propensity': '>', 'Nb_Samples': '>', 'TOP5_MUT_RANK_CI_MIXMHC': '<',
            'TOP5_MUT_RANK_CI_PRIME': '<', 'TOP5_MUT_RANK_CI_netMHCpan': '<',  'TOP5_MUT_RANK_CII': '<',
            'TOP5_WT_RANK_CII': '<', 'MIN_MUT_RANK_CII': '<', 'WT_BEST_RANK_CII': '<', 'COUNT_MUT_RANK_CII': '>',
            'rank_in_mutation': '<', 'peptide_score': '>', 'mutation_score': '>', 'mutation_rank': '<',
            'DAI_NetMHC': '<', 'DAI_MixMHC': '<', 'DAI_NetStab': '<', 'DAI_MixMHC_mbp': '<', 'DAI_aa_coeff': '>'
        }
        self.cat_features = \
            ['aa_mutant', 'Clonality', 'Zygosity', #'INTRACELLULAR_LOCATIONS',
             'FATHMM.prediction..COSMIC.', 'gene_driver_Intogen', 'mutation_driver_statement_Intogen',
             'bestWTMatchType_I', 'mut_is_binding_pos']

        self.ordinal_features = \
            ['rnaseq_gene_expression_quartile', 'mutant_other_significant_alleles', 'number_overlaping_HLA_II',
             'pep_mut_end', 'rank_in_mutation', 'number_included_HLA_I', 'mutation_rank', 'pep_mut_start',
             'pep_mut_start_9', 'pep_mut_start_10', 'pep_mut_start_11', 'pep_mut_start_12', 'seq_len']

    def get_ipmsdb_file(self):
        return self.ipmsdb_file

    def get_human_fasta_file(self):
        return self.human_fasta_file

    def get_virus_fasta_file(self):
        return self.virus_fasta_file

    def get_human_allele_score_dist_file(self):
        return self.human_allele_score_dist_file

    def get_virus_allele_score_dist_file(self):
        return self.virus_allele_score_dist_file

    def get_allotype_file(self):
        return self.allo_file

    def get_result_dir(self):
        return self.result_dir

    def get_plot_dir(self):
        return self.plot_dir

    def get_data_dir(self):
        return self.data_dir

    def get_pickle_dir(self):
        return self.pickle_dir

    def get_exe_dir(self):
        return self.exe_dir

    def patient_exist(self, patient):
        if 'Rosenberg_' in patient:
            patient = 'Rosenberg'
        return patient in self.patients

    def get_original_file(self, patient):
        if 'Rosenberg_' in patient:
            patient = 'Rosenberg'
        return os.path.join(self.data_dir, self.org_file_dict[patient])

    def get_neodisc_file(self, patient):
        if 'Rosenberg_' in patient:
            patient = 'Rosenberg'
        return os.path.join(self.data_dir, self.neodisc_file_dict[patient])

    def get_patients_with_immo(self):
        return self.patients_with_immo

    def get_patients_without_immo(self):
        return self.patients_without_immo

    def get_features(self):
        return self.features

    def get_ml_features(self):
        return self.ml_features

    def get_numerical_features(self):
        return self.num_features

    def get_categorical_features(self):
        return self.cat_features

    def get_ordinal_features(self):
        return self.ordinal_features

    def get_order_relation(self, feature):
        assert feature in self.num_features_order
        return self.num_features_order[feature]

    def get_AA(self):
        return self.aas

    def get_htide_immuno_info_file(self):
        return self.htide_info_file

    def get_gartner_info_files(self, peptide_type='long'):
        if peptide_type == 'long':
            return self.gartner_info_long_train, self.gartner_info_long_test
        else:
            return self.gartner_info_short_train, self.gartner_info_short_test

    def get_gartner_info_excel_file(self):
        return self.gartner_info_file

    def get_parkhurst_info_file(self):
        return self.parkhurst_info_file_long

    def get_tesla_info_files(self):
        return self.tesla_info_short_1, self.tesla_info_short_2

    def get_htide_patients(self):
        return self.htide_patients

    def get_rosenberg_train_patients(self):
        return self.rosenberg_train_patients

    def get_rosenberg_test_patients(self):
        return self.rosenberg_test_patients

    def get_rosenberg_patients(self):
        return self.rosenberg_train_patients + self.rosenberg_test_patients

    def get_env_variables(self):
        return self.env

    def get_protein_seq_file(self, version='37'):
        if version == '37':
            return self.protein_seq_file_37
        else:
            return self.protein_seq_file_38

    def get_data_validity_file(self):
        return self.data_validity_file

    def get_immunogenicity_info_file(self):
        return self.data_immunogenicity_info_file

    def get_cat_to_num_info_file(self, patient_set, peptide_type='long'):
        # backwards compatibility to old dataset names
        if patient_set == 'Gartner_train' or patient_set == 'NCI_train':
            return self.cat_to_num_info_files[peptide_type]['NCI_train']
        elif patient_set == 'Gartner' or patient_set == 'NCI':
            return self.cat_to_num_info_files[peptide_type]['NCI']
        else:
            return self.cat_to_num_info_files[peptide_type]['NCI']

    def get_gartner_long_result_file(self):
        return self.gartner_long_results

    def get_tesla_result_file(self):
        return self.tesla_results
