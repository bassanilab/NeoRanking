from DataWrangling.RosenbergImmunogenicityAnnotatorLong import *
from DataWrangling.RosenbergImmunogenicityAnnotatorShort import *
from Utils.DataManager import *

rosenberg_info_file = Parameters().get_gartner_info_excel_file()
data_test = pd.read_excel(open(rosenberg_info_file, 'rb'), sheet_name='Supplementary Table 18', header=1, nrows=26)
data_train = pd.read_excel(open(rosenberg_info_file, 'rb'), sheet_name='Supplementary Table 1', header=1, nrows=70)

gartner_data = pd.concat([data_test, data_train], ignore_index=True)
gartner_data = gartner_data[['ID', 'Total nmers', 'nmers screened', 'CD8+ nmers', 'Total mmps', 'Total mmps Screened',
                             'Total Confirmed mmp + Restriction Element confirmed', 'Patient HLAs']]
gartner_data.rename({'ID': 'Patient', 'Total nmers': 'Gartner et al. mutations',
                     'nmers screened': 'Gartner et al. mutations screened with minigenes',
                     'CD8+ nmers': 'Gartner et al. CD8+ immunogenic mutations screened with minigenes',
                     'Total mmps': 'Gartner et al. 8-12 mers covering mutation',
                     'Total mmps Screened': 'Gartner et al. 8-12 mers covering mutation screened with minigenes',
                     'Total Confirmed mmp + Restriction Element confirmed':
                         'Gartner et al. CD8+ immunogenic 8-12 mers covering mutation screened with minigenes',
                     'Patient HLAs': 'Gartner et al. class I alleles'}, axis=1, inplace=True)

gartner_data = gartner_data.astype({'Patient': int,
                                    'Gartner et al. mutations': int,
                                    'Gartner et al. mutations screened with minigenes': int,
                                    'Gartner et al. CD8+ immunogenic mutations screened with minigenes': int,
                                    'Gartner et al. 8-12 mers covering mutation': int,
                                    'Gartner et al. 8-12 mers covering mutation screened with minigenes': int,
                                    'Gartner et al. CD8+ immunogenic 8-12 mers covering mutation screened with minigenes': int,
                                    'Gartner et al. class I alleles': str})

mgr = DataManager(rna_seq=False, netmhc=False)
annotator_long = RosenbergImmunogenicityAnnotatorLong(mgr)
annotator_short = RosenbergImmunogenicityAnnotatorShort(mgr)

neodisc_info = []
for p in gartner_data['Patient']:
    p = str(p)
    alleles = mgr.get_classI_allotypes(p)
    alleles = list(map(lambda a: 'HLA-'+a.replace('*', ''), alleles))

    if p in mgr.get_valid_patients(peptide_type='long') and p in \
            set.union(annotator_long.get_patients('gartner_test'), annotator_long.get_patients('gartner_train')) and \
            p in mgr.get_valid_patients(peptide_type='short') and p in \
            set.union(annotator_short.get_patients('gartner_test'), annotator_short.get_patients('gartner_train')):

        neodisc_data_long = annotator_long.annotate_gartner(p)
        CD8_pos_cnt_long = sum(neodisc_data_long['response_type'] == 'CD8')
        screened_cnt_long = \
            sum((neodisc_data_long['response_type'] == 'CD8') | (neodisc_data_long['response_type'] == 'negative'))
        has_rna_seq = 'rnaseq_TPM' in neodisc_data_long.columns

        neodisc_data_short = annotator_short.annotate_patient(p)
        CD8_pos_cnt_short = sum(neodisc_data_short['response_type'] == 'CD8')
        screened_cnt_short = \
            sum((neodisc_data_short['response_type'] == 'CD8') | (neodisc_data_short['response_type'] == 'negative'))

        neodisc_info.append([p, neodisc_data_long.shape[0], screened_cnt_long, CD8_pos_cnt_long,
                             neodisc_data_short.shape[0], screened_cnt_short, CD8_pos_cnt_short, ",".join(alleles),
                             has_rna_seq])


neodisc_data = \
    pd.DataFrame(neodisc_info, columns=['Patient',
                                        'NeoDisc mutations',
                                        'NeoDisc mutations screened with minigenes',
                                        'NeoDisc CD8+ immunogenic mutations screened with minigenes',
                                        'NeoDisc 8-12 mers covering mutation',
                                        'NeoDisc 8-12 mers covering mutation screened with minigenes',
                                        'NeoDisc CD8+ immunogenic 8-12 mers covering mutation screened with minigenes',
                                        'NeoDisc class I alleles',
                                        'RNA-seq data available on dbgap'])

neodisc_data = neodisc_data.astype({'Patient': int,
                                    'NeoDisc mutations': int,
                                    'NeoDisc mutations screened with minigenes': int,
                                    'NeoDisc CD8+ immunogenic mutations screened with minigenes': int,
                                    'NeoDisc 8-12 mers covering mutation': int,
                                    'NeoDisc 8-12 mers covering mutation screened with minigenes': int,
                                    'NeoDisc CD8+ immunogenic 8-12 mers covering mutation screened with minigenes': int,
                                    'NeoDisc class I alleles': str,
                                    'RNA-seq data available on dbgap': bool})


comb_data = pd.merge(gartner_data, neodisc_data, how="inner", on=["Patient"])


def count_alleles(allele_str):
    return len(allele_str.split(','))


def count_ovrlp(allele_str1, allele_str2):
    als = allele_str1.split(',')
    cnt = 0
    for a in als:
        if a in allele_str2:
            cnt += 1
    return cnt


comb_data['NeoDisc class I allele count'] = \
    comb_data.apply(lambda r: len(r['NeoDisc class I alleles'].split(',')), axis=1)
comb_data['Gartner et al. class I allele count'] = \
    comb_data.apply(lambda r: len(r['Gartner et al. class I alleles'].split(',')), axis=1)
comb_data['NeoDisc Gartner allele ovrlp'] = \
    comb_data.apply(lambda r: count_ovrlp(r['NeoDisc class I alleles'], r['Gartner et al. class I alleles']), axis=1)

ordered_columns = \
    ['Patient',
     'Gartner et al. mutations', 'NeoDisc mutations',
     'Gartner et al. mutations screened with minigenes', 'NeoDisc mutations screened with minigenes',
     'Gartner et al. CD8+ immunogenic mutations screened with minigenes', 'NeoDisc CD8+ immunogenic mutations screened with minigenes',
     'Gartner et al. 8-12 mers covering mutation', 'NeoDisc 8-12 mers covering mutation',
     'Gartner et al. 8-12 mers covering mutation screened with minigenes',
     'NeoDisc 8-12 mers covering mutation screened with minigenes',
     'Gartner et al. CD8+ immunogenic 8-12 mers covering mutation screened with minigenes',
     'NeoDisc CD8+ immunogenic 8-12 mers covering mutation screened with minigenes',
     'Gartner et al. class I alleles', 'NeoDisc class I alleles',
     'Gartner et al. class I allele count', 'NeoDisc class I allele count', 'NeoDisc Gartner allele ovrlp',
     'RNA-seq data available on dbgap']

comb_data = comb_data[ordered_columns]

comb_data.to_csv(os.path.join(Parameters().get_plot_dir(), 'Compare_NeoDisc_Gartner.txt'),
                 header=True, index=False, sep="\t")




