import os.path

import scipy.stats as stats
from Bio import SeqIO
import math
from DataWrangling.DataLoader import *
from collections import Counter
from Features.BindingInfo.SequenceLogo import *
from sklearn.linear_model import LogisticRegression, LinearRegression


class CalcBindingAffinityScoreDists:

    def __init__(self, seq_logo_dir, seq_fasta_file, proteins=set()):

        self.sequenceLogoMgr = SequenceLogoMgr(seq_logo_dir)
        file = open(seq_fasta_file)
        if len(proteins) > 0:
            fasta_sequences = (r for r in SeqIO.parse(file, "fasta") if r.id in proteins)
        else:
            fasta_sequences = SeqIO.parse(file, 'fasta')

        self.fasta_sequences = {}
        for fasta_item in fasta_sequences:
            name, sequence = fasta_item.id, str(fasta_item.seq)
            self.fasta_sequences[name] = sequence

        file.close()

    def match_fasta(self, allele, length):

        sequenceLogo = self.sequenceLogoMgr.get_sequence_logo(allele)
        data = sequenceLogo.get_freq_data(length)

        scores = np.array([])
        for p in self.fasta_sequences.keys():
            scores = np.append(scores,
                               CalcBindingAffinityScoreDists.match_sequence(data, length, self.fasta_sequences[p]))

        return scores

    @staticmethod
    def match_sequence(data, length, protein_sequence):

        scores = np.zeros(len(protein_sequence) - length + 1)

        for i in np.arange(start=0, stop=len(protein_sequence) - length + 1):
            sum_score = 0
            for j in np.arange(length):
                aa = protein_sequence[i + j]
                if aa not in ['X', 'U', 'J', '-', '*']:
                    sum_score += data.loc[aa, j]

            scores[i] = sum_score

        return scores

    def get_number_aa(self):
        cnt = 0
        for p in self.fasta_sequences.keys():
            cnt += len(self.fasta_sequences[p])

        return cnt

    @staticmethod
    def return_block_maxima(x, n_blocks):
        np.random.shuffle(x)
        block_size = int(np.floor(len(x) / n_blocks))

        maxima = []
        for i in np.arange(n_blocks):
            start = i * block_size
            maxima = np.append(maxima, np.max(x[start:(start + block_size)]))
        return maxima

    @staticmethod
    def fit_gevd(x, n_blocks):
        maxima = CalcBindingAffinityScoreDists.return_block_maxima(x, n_blocks)
        evd_parameters = stats.genextreme.fit(maxima)
        return evd_parameters

    @staticmethod
    def get_header():
        return "Allele\tLength\tCount\tMin\tMax\tMean\tMedian\tVariance\tSkewness\tKurtosis\tPercentile_90\t" \
               "Percentile_95\tPercentile_99\tT_sigma\tT_2sigma\ttail_index\tgevd_max_loc\tgevd_max_scale\t" \
               "tail_length\n"

    @staticmethod
    def calc_dist_stats(allele, scores, peptide_length):
        ss = stats.describe(scores)
        percentiles = stats.scoreatpercentile(scores, [50, 90, 95, 99])
        evd_params = CalcBindingAffinityScoreDists.fit_gevd(scores, 1000)
        values = [allele, peptide_length, len(scores), ss.minmax[0], ss.minmax[1], ss.mean,
                  percentiles[0], ss.variance, ss.skewness, ss.kurtosis, percentiles[1], percentiles[2],
                  percentiles[3], stats.percentileofscore(scores, ss.mean + math.sqrt(ss.variance)),
                  stats.percentileofscore(scores, ss.mean + 2 * math.sqrt(ss.variance)), evd_params[0], evd_params[1],
                  evd_params[2], ss.minmax[1] - percentiles[3]]

        return values


class CalcAllelePropensity:

    def __init__(self, score_dist_stats_file, patients, mutation_types, immunogenic_types, response_types,
                 input_file_tag, peptide_type='long'):

        self.patients = patients
        self.score_dist_stats_file = score_dist_stats_file

        self.immunogenicity_ratios = \
            self.calc_immunogenicity_ratios(mutation_types, response_types, immunogenic_types, input_file_tag,
                                            peptide_type)
        self.promiscuity = self.calc_promiscuity()
        self.propensities = self.calc_propensities()

    def calc_promiscuity(self):
        promiscuity_paper_data = \
            pd.read_csv(os.path.join(Parameters().get_data_dir(), 'hla', 'promiscuity.txt'), sep="\t", header=0)
        score_dist_data = pd.read_csv(self.score_dist_stats_file, sep="\t", header=0)
        prom_alleles = np.unique(promiscuity_paper_data['Allele'])

        promiscuities = {}
        for peptide_length in [9, 10, 11, 12]:
            data_len = score_dist_data.loc[score_dist_data['Length'] == peptide_length, :]
            r = []
            X = np.empty((0, len(score_dist_data.columns) - 3), int)
            alleles = []
            for i in data_len.index:
                a = str(data_len.loc[i, 'Allele'])
                a = a.replace("*", "").replace(":", "")

                if a in prom_alleles:
                    r = np.append(r, promiscuity_paper_data.loc[promiscuity_paper_data['Allele'] == a, 'Promiscuity'])
                    X = np.append(X, [np.array(data_len.loc[i, :])[3:]], axis=0)
                    alleles = np.append(alleles, a)

            reg = LinearRegression().fit(X, r)

            X = np.empty((0, len(score_dist_data.columns) - 3), int)
            alleles = []
            for i in data_len.index:
                a = 'HLA-' + str(data_len.loc[i, 'Allele'])
                alleles = np.append(alleles, a)
                X = np.append(X, [np.array(data_len.loc[i, :])[3:]], axis=0)

            p_pred = reg.predict(X)

            pred_dict = {}
            for i in np.arange(len(alleles)):
                pred_dict[alleles[i]] = p_pred[i]

            promiscuities[peptide_length] = pred_dict

        return promiscuities

    def calc_immunogenicity_ratios(self, mutation_types, response_types, immunogenic_types, input_file_tag, peptide_type):

        features = ['response_type', 'mut_allele_0']

        data_loader = DataLoader(transformer=None, normalizer=None, features=features,
                                 mutation_types=mutation_types, response_types=response_types,
                                 immunogenic=immunogenic_types, min_nr_immuno=0)

        data, X, y = data_loader.load_patients(self.patients, input_file_tag, peptide_type)

        if data is None or data.shape[0] == 0:
            return {}

        pos_alleles = []
        neg_alleles = []
        for i in data.index:
            if data.loc[i, 'response_type'] in immunogenic_types:
                pos_alleles = np.append(pos_alleles, data.loc[i, 'mut_allele_0'])
            else:
                neg_alleles = np.append(neg_alleles, data.loc[i, 'mut_allele_0'])

        pos_cnt = Counter(pos_alleles)
        neg_cnt = Counter(neg_alleles)
        tot_cnt = pos_cnt + neg_cnt
        ratios = {}
        for (allele, cnt) in tot_cnt.items():
            ratios[allele] = pos_cnt[allele] / cnt

        return dict(sorted(ratios.items(), key=lambda item: item[1], reverse=True))

    def calc_propensities(self):

        if len(self.immunogenicity_ratios) == 0:
            return {}

        score_dist_data = pd.read_csv(self.score_dist_stats_file, sep="\t", header=0)

        propensities = {}
        for peptide_length in [9, 10, 11, 12]:
            data_len = score_dist_data.loc[score_dist_data['Length'] == peptide_length, :]
            r = []
            X = np.empty((0, len(score_dist_data.columns) - 3), int)
            alleles = []
            for i in data_len.index:
                a = 'HLA-' + str(data_len.loc[i, 'Allele'])

                if a in self.immunogenicity_ratios:
                    r = np.append(r, self.immunogenicity_ratios[a])
                    X = np.append(X, [np.array(data_len.loc[i, :])[3:]], axis=0)
                    alleles = np.append(alleles, a)

            reg = LogisticRegression(penalty='l2', C=1, class_weight='balanced'). \
                fit(X, [1 if rr > 0 else 0 for rr in r])

            X = np.empty((0, len(score_dist_data.columns) - 3), int)
            alleles = []
            for i in data_len.index:
                a = 'HLA-' + str(data_len.loc[i, 'Allele'])
                alleles = np.append(alleles, a)
                X = np.append(X, [np.array(data_len.loc[i, :])[3:]], axis=0)

            p_pred = reg.predict_proba(X)[:, 1]

            pred_dict = {}
            for i in np.arange(len(alleles)):
                pred_dict[alleles[i]] = p_pred[i]

            propensities[peptide_length] = pred_dict

        return propensities

    def process_alleles_long(self, mut_alleles, wt_alleles, peptide_lengths):

        mut_allele_propensities = []
        wt_allele_propensities = []
        for i in np.arange(len(mut_alleles)):
            mut_allele_propensities.append(self.get_propensity(mut_alleles[i], peptide_lengths[i]))
            wt_allele_propensities.append(self.get_propensity(wt_alleles[i], peptide_lengths[i]))

        return mut_allele_propensities, wt_allele_propensities

    def process_alleles_short(self, mut_alleles, peptide_lengths):
        mut_allele_propensities = []
        for i in range(len(mut_alleles)):
            allele = mut_alleles[i]
            pept_len = peptide_lengths[i]
            if ',' in allele:
                allele_list = str(allele).split(',')
                best_score = -1000000
                for a in allele_list:
                    af = CalcAllelePropensity.format_allele_short(a)
                    score = self.get_propensity(af, pept_len)
                    if score > best_score:
                        best_score = score

                mut_allele_propensities.append(best_score)
            else:
                allele = CalcAllelePropensity.format_allele_short(allele)
                mut_allele_propensities.append(self.get_propensity(allele, pept_len))

        return mut_allele_propensities

    def get_propensity(self, allele, peptide_length):
        if peptide_length <= 8:
            peptide_length = 9
        if peptide_length > 12:
            peptide_length = 12

        if allele in self.propensities[peptide_length]:
            return self.propensities[peptide_length][allele]
        else:
            return 0.0

    def get_immunogenicity_ratios(self):
        return self.immunogenicity_ratios

    @staticmethod
    def merge_with_data_long(data, mut_allele_propensities, wt_allele_propensities, i):
        data['mut_allele_propensity_' + str(i)] = mut_allele_propensities
        data['wt_allele_propensity_' + str(i)] = wt_allele_propensities
        return data

    @staticmethod
    def merge_with_data_short(data, mut_allele_propensities):
        data['mut_allele_propensity'] = mut_allele_propensities
        return data

    def add_features(self, patient, file_tag_input, file_tag_output, peptide_type='long', write_res=True):
        mgr = DataManager()
        parameters = Parameters()

        data = mgr.get_processed_data(patient, file_tag_input, peptide_type)

        if peptide_type == 'long':
            if 'mut_peptide_0' not in data.columns:
                print('Patient '+patient+' does not contain netmhc predictions. Unable to add binding info.')
                return None
            data = self.add_features_long(data)
        else:
            if 'mutant_best_alleles_netMHCpan' not in data.columns:
                print('Patient '+patient+' does not contain netmhc prediction. Unable to add binding info.')
                return None
            data = self.add_features_short(data)

        mgr.put_processed_data(data, patient, file_tag_output, peptide_type)

        if write_res:
            out_file = os.path.join(parameters.get_result_dir(), patient+"_"+peptide_type+"_"+file_tag_output+".txt")
            data.to_csv(path_or_buf=out_file, sep="\t", index=False, header=True)

        return data

    def add_features_long(self, data):
        netMHCpan_ranks = [int(c[c.rfind('_') + 1:]) for c in data.columns if 'mut_peptide_pos_' in c]
        for i in netMHCpan_ranks:
            mut_alleles, wt_alleles, peptide_lengths = CalcAllelePropensity.get_peptide_info_long(data, i)
            mut_alleles_prop, wt_alleles_prop = self.process_alleles_long(mut_alleles, wt_alleles, peptide_lengths)
            data = CalcAllelePropensity.merge_with_data_long(data, mut_alleles_prop, wt_alleles_prop, i)

        return data

    def add_features_short(self, data):
        mut_alleles, peptide_lengths = CalcAllelePropensity.get_peptide_info_short(data)
        mut_alleles_prop = self.process_alleles_short(mut_alleles, peptide_lengths)
        data = CalcAllelePropensity.merge_with_data_short(data, mut_alleles_prop)

        return data

    @staticmethod
    def get_peptide_info_long(data, index):
        mut_allele_column = 'mut_allele_' + str(index)
        wt_allele_column = 'wt_allele_' + str(index)
        mut_peptide_column = 'mut_peptide_' + str(index)
        if mut_allele_column not in data.columns:
            return None, None, None
        else:
            peptide_lengths = [len(str(p)) for p in data[mut_peptide_column]]
            return np.array(data[mut_allele_column], dtype='str'), np.array(data[wt_allele_column], dtype='str'), \
                   np.array(peptide_lengths, dtype='int')

    @staticmethod
    def get_peptide_info_short(data):
        mut_allele_column = 'mutant_best_alleles_netMHCpan'
        mut_peptide_column = 'mutant_seq'
        if mut_allele_column not in data.columns:
            return None, None
        else:
            peptide_lengths = [len(str(p)) for p in data[mut_peptide_column]]
            return np.array(data[mut_allele_column], dtype='str'), np.array(peptide_lengths, dtype='int')

    @staticmethod
    def format_allele_short(allele):
        return 'HLA-'+allele[0]+'*'+allele[1:3]+':'+allele[3:]
