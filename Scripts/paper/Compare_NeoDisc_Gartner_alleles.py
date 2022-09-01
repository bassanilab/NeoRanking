from DataWrangling.RosenbergImmunogenicityAnnotatorLong import *
from DataWrangling.RosenbergImmunogenicityAnnotatorShort import *
from Utils.DataManager import *

rosenberg_info_file = Parameters().get_gartner_info_excel_file()
data_test = pd.read_excel(open(rosenberg_info_file, 'rb'), sheet_name='Supplementary Table 18', header=1, nrows=26)
data_train = pd.read_excel(open(rosenberg_info_file, 'rb'), sheet_name='Supplementary Table 1', header=1, nrows=70)

gartner_data = pd.concat([data_test, data_train], ignore_index=True)
gartner_data = gartner_data[['ID', 'Patient HLAs']]

gartner_data.rename({'ID': 'Patient', 'Patient HLAs': 'Gartner class I alleles'}, axis=1, inplace=True)
gartner_data = gartner_data.astype({'Patient': str, 'Gartner class I alleles': str})

mgr = DataManager()

neodisc_info = []
for p in gartner_data['Patient']:
    if p in mgr.get_valid_patients(peptide_type='long'):
        alleles = mgr.get_classI_allotypes(str(p))
        if len(alleles) > 0:
            alleles = list(map(lambda a: 'HLA-'+a.replace('*', ''), alleles))
            neodisc_info.append([p, ", ".join(alleles)])

neodisc_data = pd.DataFrame(neodisc_info, columns=['Patient', 'NeoDisc class I alleles'])
neodisc_data = neodisc_data.astype({'Patient': str, 'NeoDisc class I alleles': str})

comb_data = pd.merge(gartner_data, neodisc_data, how="inner", on=["Patient"])


def count_alleles(allele_str):
    return len(allele_str.split(', '))


def count_allele_overlap(allele_str1, allele_str2):
    als = allele_str1.split(', ')
    cnt = 0
    for a in als:
        if a in allele_str2:
            cnt += 1
    return cnt


comb_data['NeoDisc class I allele count'] = \
    comb_data.apply(lambda r: len(r['NeoDisc class I alleles'].split(',')), axis=1)
comb_data['Gartner class I allele count'] = \
    comb_data.apply(lambda r: len(r['Gartner class I alleles'].split(',')), axis=1)
comb_data['NeoDisc-Gartner allele overlap count'] = \
    comb_data.apply(lambda r: count_allele_overlap(r['NeoDisc class I alleles'], r['Gartner class I alleles']), axis=1)

ordered_columns = \
    ['Patient',
     'Gartner class I alleles',
     'NeoDisc class I alleles',
     'Gartner class I allele count', 
     'NeoDisc class I allele count', 
     'NeoDisc-Gartner allele overlap count']

comb_data = comb_data[ordered_columns]

comb_data.to_csv(os.path.join(Parameters().get_plot_dir(), 'Compare_NeoDisc_Gartner_alleles.txt'),
                 header=True, index=False, sep="\t")

cnt = comb_data.shape[0]
idx = comb_data.apply(lambda row:
                      row['NeoDisc class I allele count'] == row['NeoDisc-Gartner allele overlap count'] and
                      row['NeoDisc class I allele count'] == row['Gartner class I allele count'],
                      axis=1)
cnt_match = sum(idx)

print("From {0} patients, {1} have same alleles".format(cnt, cnt_match))
no_match = comb_data.loc[~ np.array(idx), :]
no_match.apply(lambda row:
               print("Patient: {0}\tGartner: {1}\tNeoDisc: {2}".
                     format(row['Patient'],
                            set.difference(set(row['Gartner class I alleles'].split(', ')),
                                           set(row['NeoDisc class I alleles'].split(', '))),
                            set.difference(set(row['NeoDisc class I alleles'].split(', ')),
                                           set(row['Gartner class I alleles'].split(', '))))), axis=1)



