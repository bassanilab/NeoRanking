import os
import os.path
from typing import Final


class GlobalParameters:
    """
    Class to store global parameters
    Attributes:
        base_dir (str): base directory of project files
        data_dir (str): directory that holds data files
        plot_dir (str): directory that holds figure files
        classifier_result_dir (str): directory that holds classifier result files
        classifier_model_dir (str): directory that holds classifier model files
        neopep_data_org_file (str): tab file containing all neo-peptide data
        mutation_data_org_file (str): tab file containing all mutation data
        neopep_data_ml_sel_file (str): tab file containing rows of neo-peptide data selected for ML
        mutation_data_ml_sel_file (str): tab file containing rows of mutation data selected for ML
        neopep_data_ml_file (str): tab file containing neo-peptide data normalized for ML
        mutation_data_ml_file (str): tab file containing mutation data data normalized for ML
        neopep_data_plot_file (str): tab file containing neo-peptide data normalized for histogram and scatter plots
        mutation_data_plot_file (str): tab file containing mutation data normalized for histogram and scatter plots
        cat_to_num_info_files (dict[str, dict[str, str]]): dictionary with file names for imputation of categorical
                                                           variables
        tesla_result_file (str): results from TESLA paper containing FR, TTIF, and AUPRC scores of different groups
        gartner_nmer_train_file (str): training data matrix from Gartner et al with mutation features and immunogenicity
                                       annotation downloaded from figshare link provided in Gartner et al
        gartner_nmer_test_file (str): testing data matrix from Gartner et al with mutation features and immunogenicity
                                       annotation downloaded from figshare link provided in Gartner et al
        gartner_nmer_rank_file (str): file containing the ranking of mutations in NCI_test obtained by Gartner et al.
        gartner_mmp_rank_file (str): file containing the ranking of neo-peptides in NCI_test obtained by Gartner et al.
        hlaI_allele_file (str): file containing the HLA class I alleles of all patients
        datasets (list[str]): datasets used in this study  ['NCI', 'NCI_train', 'NCI_test', 'TESLA', 'HiTIDE']
        datasets_encoding (list[str]): datasets used for encoding categorical values  ['NCI', 'NCI_train']
        peptide_types (list[str]): peptide types ['neopep', 'mutation']
        objectives (list[str]): objectives for data normalization ['ml', 'plot']
        response_types (list[str]): immunogenicity measurement response types ['CD8', 'negative', 'not_tested']
        mutation_types (list[str]): mutation types to include ['SNV', 'INSERTION', 'DELETION', 'FSS']
        classifiers (list[str]): classifiers used in this study
        aas (list[str]): list of amino acids
        ml_features_neopep (list[str]): list of features used for classification of neo-peptides
        features_neopep (list[str]): list of features for neo-peptides
        feature_types_neopep (dict[str, any]): types of features_neopep
        ml_feature_mv_neopep (dict[str, str]): order of features_neopep values (used for missing value imputation)
        ml_features_mutation (list[str]): list of features used for classification of neo-peptides
        features_mutation (list[str]): list of features for neo-peptides
        feature_types_mutation (dict[str, any]): types of features_neopep
        ml_feature_mv_mutation (dict[str, str]): order of features_mutation values (used for missing value imputation)
        nr_hyperopt_rep (int): number of replicate hyperopt runs
        nr_hyperopt_iter (int): number of hyperopt iterations
        nr_hyperopt_cv (int): number of hyperopt cross-validation folds
        neopep_alpha (float): value of alpha in rank_score function used for training neo-peptides
        mutation_alpha (float): value of alpha in rank_score function used for training mutations
        normalizer (str): normalizer to be used ('q': quantile, 'p': power, 'z': standard, 'i': minmax, 'l': log, 'a': asinh, 'n': none)
        nr_non_immuno_neopeps (int): nr non-immunogenic peptides sampled
        cat_type (str): conversion of categorical to numerical values. either 'float' or 'int'
        max_netmhc_rank (float): maximal netmhc rank for neo-peptide. -1 if no filter applied
        excluded_genes (list): peptides of these genes are excluded from prioritization
        plot_normalization (dict): feature normalization for plots only (not for ML)
        plot_feature_names (dict): feature names used in plots
        color_immunogenic (str): color used to represent immunogenic peptides in plots
        color_negative (str): color used to represent non-immunogenic peptides in plots
    """

    base_dir: Final[str] = os.getenv('NEORANKING_RESOURCE')
    data_dir: Final[str] = os.path.join(base_dir, "data")
    plot_dir: Final[str] = os.path.join(base_dir, "plots")
    classifier_result_dir: Final[str] = os.path.join(base_dir, "classifier_results")
    classifier_model_dir: Final[str] = os.path.join(base_dir, "classifier_models")

    neopep_data_org_file: Final[str] = os.path.join(data_dir, "Neopep_data_org.txt")
    mutation_data_org_file: Final[str] = os.path.join(data_dir, "Mutation_data_org.txt")
    neopep_data_ml_sel_file: Final[str] = os.path.join(data_dir, "Neopep_data_ml_sel.txt")
    mutation_data_ml_sel_file: Final[str] = os.path.join(data_dir, "Mutation_data_ml_sel.txt")
    neopep_data_ml_file: Final[str] = os.path.join(data_dir, "Neopep_data_ml_norm.txt")
    mutation_data_ml_file: Final[str] = os.path.join(data_dir, "Mutation_data_ml_norm.txt")
    neopep_data_plot_file: Final[str] = os.path.join(data_dir, "Neopep_data_plot_norm.txt")
    mutation_data_plot_file: Final[str] = os.path.join(data_dir, "Mutation_data_plot_norm.txt")

    cat_to_num_info_files: Final[dict] = \
        {
            'neopep': {'NCI_train': os.path.join(data_dir, 'cat_encoding', 'Cat_to_num_info_neopep_NCI_train.txt'),
                       'NCI': os.path.join(data_dir, 'cat_encoding', 'Cat_to_num_info_neopep_NCI_all.txt')},
            'mutation': {'NCI_train': os.path.join(data_dir, 'cat_encoding', 'Cat_to_num_info_mutation_NCI_train.txt'),
                         'NCI': os.path.join(data_dir, 'cat_encoding', 'Cat_to_num_info_mutation_NCI_all.txt')}
        }

    tesla_result_file: Final[str] = os.path.join(data_dir, "mmc5.xlsx")
    gartner_nmer_train_file: Final[str] = os.path.join(data_dir, 'NmersTrainingSet.txt')
    gartner_nmer_test_file: Final[str] = os.path.join(data_dir, 'NmersTestingSet.txt')
    gartner_nmer_rank_file: Final[str] = os.path.join(data_dir, 'Gartner_nmers_ranking.txt')
    gartner_mmp_rank_file: Final[str] = os.path.join(data_dir, 'Gartner_mmps_ranking.txt')
    hlaI_allele_file: Final[str] = os.path.join(data_dir, 'hla', 'HLA_allotypes.txt')

    datasets: Final[list] = ['NCI', 'NCI_train', 'NCI_test', 'TESLA', 'HiTIDE']
    datasets_encoding: Final[list] = ['NCI', 'NCI_train']
    peptide_types: Final[list] = ['neopep', 'mutation']
    objectives: Final[list] = ['ml', 'plot']
    response_types: Final[list] = ['CD8', 'negative', 'not_tested']
    mutation_types: Final[list] = ['SNV', 'INSERTION', 'DELETION', 'FSS']

    aas: Final[list] = \
        ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    classifiers = ['SVM', 'SVM-lin', 'RF', 'CART', 'ADA', 'LR', 'NNN', 'XGBoost']
    neopep_alpha: Final[float] = 0.005
    mutation_alpha: Final[float] = 0.05
    nr_hyperopt_rep = 10
    nr_hyperopt_iter = 200
    nr_hyperopt_cv = 5
    normalizer: Final[str] = 'q'
    nr_non_immuno_neopeps: Final[int] = 200000
    cat_type: Final[str] = 'float'  # either float or int
    max_netmhc_rank: Final[int] = -1

    excluded_genes: Final[list] = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5',
                                   'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DMA', 'TRBV3', 'TRBV5',
                                   'TRBV6', 'TRBV6-1', 'TRBV10', 'TRBV10-1', 'TRBV11', 'TRAV12', 'KRT1', 'PRSS3']

    ml_features_neopep: Final[list] = \
        ['CCF', 'Clonality', 'rnaseq_TPM', 'rnaseq_alt_support', 'CSCAPE_score',
         'mutant_other_significant_alleles', 'mutant_rank', 'mutant_rank_PRIME',
         'mutant_rank_netMHCpan', 'Sample_Tissue_expression_GTEx',
         'GTEx_all_tissues_expression_mean', 'TCGA_Cancer_expression',
         'gene_driver_Intogen', 'nb_same_mutation_Intogen',
         'mutation_driver_statement_Intogen', 'bestWTMatchScore_I',
         'bestWTMatchOverlap_I', 'bestMutationScore_I', 'bestWTMatchType_I',
         'bestWTPeptideCount_I', 'mut_Rank_Stab', 'mut_netchop_score_ct',
         'TAP_score', 'mut_is_binding_pos', 'mut_binding_score', 'mut_aa_coeff',
         'seq_len', 'DAI_NetMHC', 'DAI_MixMHC', 'DAI_NetStab', 'DAI_MixMHC_mbp']

    features_neopep: Final[list] = \
        ['patient', 'dataset', 'train_test', 'response_type', 'Nb_Samples', 'Sample_Tissue', 'Cancer_Type',
         'chromosome', 'genomic_coord', 'ref', 'alt', 'gene', 'protein_coord', 'aa_mutant', 'aa_wt',
         'pep_mut_start', 'TumorContent', 'Zygosity', 'mutation_type'] + ml_features_neopep

    feature_types_neopep: Final[dict] = {
        'patient': 'str',
        'dataset': 'category',
        'train_test': 'category',
        'response_type': 'category',
        'Nb_Samples': 'str',
        'Sample_Tissue': 'str',
        'Cancer_Type': 'str',
        'chromosome': 'str',
        'genomic_coord': 'int64',
        'ref': 'str',
        'alt': 'str',
        'gene': 'str',
        'protein_coord': 'int32',
        'aa_mutant': 'category',
        'aa_wt': 'category',
        'mutant_seq': 'str',
        'wt_seq': 'str',
        'pep_mut_start': 'int8',
        'TumorContent': 'float64',
        'CCF': 'float64',
        'Clonality': 'category',
        'Zygosity': 'category',
        'mutation_type': 'category',
        'mutant_rank': 'float64',
        'mutant_rank_netMHCpan': 'float64',
        'mutant_rank_PRIME': 'float64',
        'mut_Rank_Stab': 'float64',
        'TAP_score': 'float64',
        'mut_netchop_score_ct': 'float64',
        'mut_binding_score': 'float64',
        'mut_is_binding_pos': 'bool',
        'mut_aa_coeff': 'float64',
        'DAI_NetMHC': 'float64',
        'DAI_MixMHC': 'float64',
        'DAI_NetStab': 'float64',
        'mutant_other_significant_alleles': 'int8',
        'DAI_MixMHC_mbp': 'float64',
        'rnaseq_TPM': 'float64',
        'rnaseq_alt_support': 'float64',
        'GTEx_all_tissues_expression_mean': 'float64',
        'Sample_Tissue_expression_GTEx': 'float64',
        'TCGA_Cancer_expression': 'float64',
        'bestWTMatchScore_I': 'float64',
        'bestWTMatchOverlap_I': 'float64',
        'bestMutationScore_I': 'float64',
        'bestWTPeptideCount_I': 'int32',
        'bestWTMatchType_I': 'category',
        'CSCAPE_score': 'float64',
        'nb_same_mutation_Intogen': 'int32',
        'mutation_driver_statement_Intogen': 'category',
        'gene_driver_Intogen': 'category',
        'seq_len': 'int8'
    }

    ml_feature_mv_neopep: Final[dict] = {
        'mutant_rank': 'max',
        'mutant_rank_netMHCpan': 'max',
        'mutant_rank_PRIME': 'max',
        'mut_Rank_Stab': 'max',
        'TAP_score': 'min',
        'mut_netchop_score_ct': 'min',
        'mut_binding_score': 'min',
        'mut_is_binding_pos': 'cnt',
        'mut_aa_coeff': 'cnt',
        'DAI_NetMHC': 'max',
        'DAI_MixMHC': 'max',
        'DAI_NetStab': 'max',
        'mutant_other_significant_alleles': 'min',
        'DAI_MixMHC_mbp': 'max',
        'rnaseq_TPM': 'min',
        'rnaseq_alt_support': 'min',
        'GTEx_all_tissues_expression_mean': 'min',
        'Sample_Tissue_expression_GTEx': 'min',
        'TCGA_Cancer_expression': 'min',
        'bestWTMatchScore_I': 'min',
        'bestWTMatchOverlap_I': 'min',
        'bestMutationScore_I': 'min',
        'bestWTPeptideCount_I': 'min',
        'CSCAPE_score': 'min',
        'CCF': 0.9,
        'nb_same_mutation_Intogen': 'min',
        'seq_len': 'cnt'
    }

    ml_features_mutation: Final[list] = \
        ['CCF', 'Clonality', 'Zygosity', 'Sample_Tissue_expression_GTEx',
         'TCGA_Cancer_expression', 'rnaseq_TPM', 'rnaseq_alt_support',
         'MIN_MUT_RANK_CI_MIXMHC', 'COUNT_MUT_RANK_CI_MIXMHC',
         'WT_BEST_RANK_CI_MIXMHC', 'MIN_MUT_RANK_CI_PRIME',
         'COUNT_MUT_RANK_CI_PRIME', 'WT_BEST_RANK_CI_PRIME',
         'COUNT_MUT_RANK_CI_netMHCpan', 'CSCAPE_score', 'gene_driver_Intogen',
         'nb_mutations_in_gene_Intogen', 'nb_same_mutation_Intogen',
         'mutation_driver_statement_Intogen', 'GTEx_all_tissues_expression_mean',
         'bestWTMatchScore_I', 'bestWTMatchOverlap_I', 'bestMutationScore_I',
         'bestWTPeptideCount_I', 'mut_Rank_EL_0', 'wt_Rank_EL_0',
         'mut_Rank_EL_1', 'wt_Rank_EL_1', 'mut_Rank_EL_2', 'wt_Rank_EL_2',
         'mut_Rank_Stab_0', 'mut_Rank_Stab_1', 'mut_Rank_Stab_2',
         'mut_netchop_score', 'mut_TAP_score_0', 'next_best_BA_mut_ranks',
         'DAI_0', 'DAI_1', 'DAI_2']

    features_mutation: Final[list] = \
        ['patient', 'dataset', 'train_test', 'response_type', 'Nb_Samples', 'Sample_Tissue', 'Cancer_Type',
         'chromosome', 'genomic_coord', 'ref', 'alt', 'gene', 'protein_coord', 'aa_mutant', 'aa_wt', 'pep_mut_start',
         'TumorContent', 'mutation_type'] + ml_features_mutation

    feature_types_mutation: Final[dict] = {
        'patient': 'category',
        'dataset': 'category',
        'train_test': 'category',
        'response_type': 'category',
        'Nb_Samples': 'str',
        'Sample_Tissue': 'str',
        'Cancer_Type': 'str',
        'chromosome': 'str',
        'genomic_coord': 'int64',
        'ref': 'str',
        'alt': 'str',
        'gene': 'str',
        'protein_coord': 'int32',
        'aa_mutant': 'category',
        'aa_wt': 'category',
        'mutant_seq': 'str',
        'wt_seq': 'str',
        'pep_mut_start': 'int8',
        'TumorContent': 'float64',
        'CCF': 'float64',
        'Clonality': 'category',
        'Zygosity': 'category',
        'mutation_type': 'category',
        'nb_same_mutation_Intogen': 'float64',
        'nb_mutations_in_gene_Intogen': 'float64',
        'mutation_driver_statement_Intogen': 'category',
        'gene_driver_Intogen': 'category',
        'rnaseq_TPM': 'float64',
        'TCGA_Cancer_expression': 'float64',
        'bestMutationScore_I': 'float64',
        'bestWTPeptideCount_I': 'int32',
        'bestWTMatchScore_I': 'float64',
        'bestWTMatchOverlap_I': 'float64',
        'rnaseq_alt_support': 'float64',
        'CSCAPE_score': 'float64',
        'GTEx_all_tissues_expression_mean': 'float64',
        'Sample_Tissue_expression_GTEx': 'float64',
        'COUNT_MUT_RANK_CI_MIXMHC': 'int32',
        'COUNT_MUT_RANK_CI_PRIME': 'int32',
        'COUNT_MUT_RANK_CI_netMHCpan': 'int32',
        'MIN_MUT_RANK_CI_MIXMHC': 'float64',
        'WT_BEST_RANK_CI_MIXMHC': 'float64',
        'MIN_MUT_RANK_CI_PRIME': 'float64',
        'WT_BEST_RANK_CI_PRIME': 'float64',
        'next_best_BA_mut_ranks': 'float64',
        'mut_Rank_EL_0': 'float64',
        'mut_Rank_EL_1': 'float64',
        'mut_Rank_EL_2': 'float64',
        'wt_Rank_EL_0': 'float64',
        'wt_Rank_EL_1': 'float64',
        'wt_Rank_EL_2': 'float64',
        'mut_Rank_Stab_0': 'float64',
        'mut_Rank_Stab_1': 'float64',
        'mut_Rank_Stab_2': 'float64',
        'DAI_0': 'float64',
        'DAI_1': 'float64',
        'DAI_2': 'float64',
        'mut_TAP_score_0': 'float64',
        'mut_netchop_score': 'float64'
    }

    ml_feature_mv_mutation: Final[dict] = {
        'nb_same_mutation_Intogen': 'min',
        'nb_mutations_in_gene_Intogen': 'min',
        'rnaseq_TPM': 'min',
        'TCGA_Cancer_expression': 'min',
        'bestMutationScore_I': 'min',
        'bestWTPeptideCount_I': 'min',
        'bestWTMatchScore_I': 'min',
        'bestWTMatchOverlap_I': 'min',
        'rnaseq_alt_support': 'min',
        'CCF': 0.9,
        'CSCAPE_score': 'min',
        'GTEx_all_tissues_expression_mean': 'min',
        'Sample_Tissue_expression_GTEx': 'min',
        'COUNT_MUT_RANK_CI_MIXMHC': 'min',
        'COUNT_MUT_RANK_CI_PRIME': 'min',
        'COUNT_MUT_RANK_CI_netMHCpan': 'min',
        'MIN_MUT_RANK_CI_MIXMHC': 'max',
        'WT_BEST_RANK_CI_MIXMHC': 'max',
        'MIN_MUT_RANK_CI_PRIME': 'max',
        'WT_BEST_RANK_CI_PRIME': 'max',
        'next_best_BA_mut_ranks': 'max',
        'mut_Rank_EL_0': 'max',
        'mut_Rank_EL_1': 'max',
        'mut_Rank_EL_2': 'max',
        'wt_Rank_EL_0': 'max',
        'wt_Rank_EL_1': 'max',
        'wt_Rank_EL_2': 'max',
        'mut_Rank_Stab_0': 'max',
        'mut_Rank_Stab_1': 'max',
        'mut_Rank_Stab_2': 'max',
        'DAI_0': 'cnt',
        'DAI_1': 'cnt',
        'DAI_2': 'cnt',
        'mut_TAP_score_0': 'min',
        'mut_netchop_score': 'min'
    }

    #
    # Visualization
    #
    color_immunogenic = 'darkorange'
    color_negative = 'royalblue'
    plot_file_formats = ['pdf', 'svg', 'png']

    plot_normalization: Final[dict] = \
        {'mutant_rank_PRIME': 'l', 'wt_best_rank_PRIME': 'l', 'mutant_rank': 'l', 'wt_best_rank': 'l',
         'mutant_rank_netMHCpan': 'l', 'wt_best_rank_netMHCpan': 'l', 'mut_Rank_Stab': 'l', 'wt_Rank_Stab': 'l',
         'mut_Stab_Score': 'n', 'wt_Stab_Score': 'n', 'TAP_score': 'n', 'mut_netchop_score_ct': 'n',
         'mut_binding_score': 'n', 'mut_is_binding_pos': 'n', 'pep_mut_start': 'i', 'mut_aa_coeff': 'n', 'DAI': 'n',
         'rnaseq_TPM': 'a', 'rnaseq_alt_support': 'n', 'GTEx_all_tissues_expression_mean': 'a',
         'Sample_Tissue_expression_GTEx': 'a', 'TCGA_Cancer_expression': 'a', 'bestWTMatchScore_I': 'a',
         'bestWTMatchOverlap_I': 'n', 'bestMutationScore_I': 'a', 'bestWTPeptideCount_I': 'a', 'bestWTMatchType_I': 'n',
         'mutant_other_significant_alleles': 'n', 'CSCAPE_score': 'n', 'Clonality': 'n',
         'CCF': 'n', 'nb_same_mutation_Intogen': 'a', 'nb_mutations_in_gene_Intogen': 'a',
         'nb_mutations_same_position_Intogen': 'a', 'mutation_driver_statement_Intogen': 'n',
         'gene_driver_Intogen': 'n', 'DAI_NetMHC': 'n', 'DAI_MixMHC': 'n', 'DAI_NetStab': 'n',
         'DAI_MixMHC_mbp': 'n', 'seq_len': 'n', 'DAI_aa_coeff': 'n', 'mut_Rank_EL_0': 'l',
         'mut_Rank_EL_1': 'l', 'mut_Rank_EL_2': 'l', 'wt_Rank_EL_0': 'l', 'wt_Rank_EL_1': 'l', 'wt_Rank_EL_2': 'l',
         'mut_Rank_Stab_0': 'l', 'mut_Rank_Stab_1': 'l', 'mut_Rank_Stab_2': 'l', 'DAI_0': 'n', 'DAI_1': 'n',
         'DAI_2': 'n', 'mut_TAP_score_0': 'n', 'mut_netchop_score': 'n', 'COUNT_MUT_RANK_CI_MIXMHC': 'n',
         'COUNT_MUT_RANK_CI_PRIME': 'n', 'COUNT_MUT_RANK_CI_netMHCpan': 'n', 'mut_nr_strong_binders_0': 'n',
         'mut_nr_weak_binding_alleles_0': 'n', 'MIN_MUT_RANK_CI_MIXMHC': 'l', 'WT_BEST_RANK_CI_MIXMHC': 'l',
         'MIN_MUT_RANK_CI_PRIME': 'l', 'WT_BEST_RANK_CI_PRIME': 'l', 'next_best_BA_mut_ranks': 'l'
         }

    plot_feature_names: Final[dict] = \
        {'mutant_rank': 'MixMHCpred Rank', 'mutant_rank_netMHCpan': 'NetMHCpan Rank', 'mutant_rank_PRIME': 'PRIME Rank',
         'mut_Rank_Stab': 'NetStab Rank', 'TAP_score': 'NetTAP Score', 'mut_netchop_score_ct': 'NetChop CT Score',
         'mut_binding_score': 'MixMHCpred Score at Mutation', 'mut_is_binding_pos': 'Mutation at Anchor',
         'pep_mut_start': 'Mutation Position', 'mut_aa_coeff': 'PRIME Coeff at Mutation',
         'DAI_NetMHC': 'NetMHCpan log_Rank DAI', 'DAI_MixMHC': 'MixMHCpred log_Rank DAI',
         'DAI_NetStab': 'NetStab log_Rank DAI', 'mutant_other_significant_alleles': 'Number Binding Alleles',
         'DAI_MixMHC_mbp': 'MixMHCpred Score DAI', 'rnaseq_TPM': 'RNAseq Expression(TPM)',
         'rnaseq_alt_support': 'RNAseq Mutation Coverage',
         'GTEx_all_tissues_expression_mean': 'GTEx Mean Tissue Expression',
         'Sample_Tissue_expression_GTEx': 'GTEx Sample Tissue Expression',
         'TCGA_Cancer_expression': 'TCGA Cancer Expression',
         'bestWTMatchScore_I': 'ipMSDB Peptide Score', 'bestWTMatchOverlap_I': 'ipMSDB Peptide Overlap',
         'bestMutationScore_I': 'ipMSDB Mutation Score', 'bestWTPeptideCount_I': 'ipMSDB Peptide Count',
         'bestWTMatchType_I': 'ipMSDB Peptide Match Type', 'CSCAPE_score': 'CSCAPE Score', 'Zygosity': 'Zygosity',
         'Clonality': 'Clonality', 'CCF': 'Cancer Cell Fraction',
         'nb_same_mutation_Intogen': 'Intogen Same Mutation Count',
         'nb_mutations_in_gene_Intogen': 'Intogen Gene Mutation Count',
         'nb_mutations_same_position_Intogen': 'Intogen Mutation Same Position Count',
         'mutation_driver_statement_Intogen': 'Intogen Mutation Driver Statement',
         'gene_driver_Intogen': 'Gene Driver Intogen', 'pep_mut_start_9': 'Mutation Position Length 9',
         'pep_mut_start_10': 'Mutation Position Length 10', 'pep_mut_start_11': 'Mutation Position Length 11',
         'pep_mut_start_12': 'Mutation Position Length 12', 'seq_len': 'Peptide Length',
         'DAI_aa_coeff': 'PRIME Coefficient DAI', 'COUNT_MUT_RANK_CI_MIXMHC': 'MixMHCpred Binding Peptide Count',
         'COUNT_MUT_RANK_CI_PRIME': 'PRIME Binding Peptide Count',
         'COUNT_MUT_RANK_CI_netMHCpan': 'NetMHC Binding Peptide Count',
         'MIN_MUT_RANK_CI_MIXMHC': 'Minimal Mut MixMHCpred Rank', 'MIN_MUT_RANK_CI_PRIME': 'Minimal Mut PRIME Rank',
         'WT_BEST_RANK_CI_MIXMHC': 'Minimal WT MixMHCpred Rank', 'WT_BEST_RANK_CI_PRIME': 'Minimal WT PRIME Rank',
         'next_best_BA_mut_ranks': 'Second Mut BA rank', 'mut_Rank_EL_0': 'Best Mut EL Rank',
         'mut_Rank_EL_1': 'Second Mut EL Rank', 'mut_Rank_EL_2': 'Third Mut EL Rank', 'wt_Rank_EL_0': 'Best WT EL Rank',
         'wt_Rank_EL_1': 'Second WT EL Rank', 'wt_Rank_EL_2': 'Third WT EL Rank',
         'mut_Rank_Stab_0': 'Best Mut Stab Rank',
         'mut_Rank_Stab_1': 'Second Mut Stab Rank', 'mut_Rank_Stab_2': 'Third Mut Stab Rank',
         'DAI_0': 'BEST EL Rank DAI',
         'DAI_1': 'Second EL Rank DAI', 'DAI_2': 'Third EL Rank DAI', 'mut_TAP_score_0': 'Best Mut TAP Score',
         'mut_netchop_score': 'Best Mut NetChop Score'
         }

    @staticmethod
    def get_cat_to_num_info_file(dataset: str, peptide_type: str):
        if dataset in GlobalParameters.datasets_encoding:
            return GlobalParameters.cat_to_num_info_files[peptide_type][dataset]
        else:
            return None
