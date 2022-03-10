import tempfile
import subprocess
import multiprocessing
from abc import abstractmethod
import os
from Utils.Parameters import *

import numpy as np


class NetMHCRunner:

    def __init__(self, score_key='Rank_EL', operator='<'):
        self.score_key = score_key
        self.operator = operator  # '<' : the smaller the better, '>' : the larger the better

        return

    @abstractmethod
    def build_cmd(self, peptide_file, peptide_file_type, allele, length, show_cmd=False):
        return ""

    def process_peptides(self, indexes, peptides, input_type, alleles, length=9, show_cmd=True):
        """ Ĉreate fasta file for with peptides and run netMHCpan """

        if input_type == 'fa':
            pep_file = tempfile.mktemp() + '.fa'
            with open(pep_file, 'w') as f:
                for i in indexes:
                    if len(peptides[i]) > 7:
                        f.write('>' + str(i) + '\n' + peptides[i] + '\n')
        else:
            pep_file = tempfile.mktemp() + '.txt'
            with open(pep_file, 'w') as f:
                for i in indexes:
                    if len(peptides[i]) > 7:
                        f.write(peptides[i] + '\n')

        with multiprocessing.Pool(len(alleles)) as pool:
            result_async = [pool.apply_async(self.process_allele, args=(pep_file, input_type, a, length, show_cmd))
                            for a in alleles]
            res_list = [r.get() for r in result_async]

        os.remove(pep_file)

        return res_list

    def format_allele(self, allele):
        return allele.replace('*', '')

    def process_allele(self, peptide_file, peptide_file_type, allele, length=9, show_cmd=True):
        """Define and call netMHCpan command line."""

        allele = self.format_allele(allele)

        cmd = self.build_cmd(peptide_file, peptide_file_type, allele, length, show_cmd)

        try:
            my_env = Parameters().get_env_variables()
            temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash', env=my_env)
#            temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
        except Exception as e:
            print(e)
            return

        res = self.read_result(temp)

        if res is None or res.shape[0] == 0:
            print('empty result, check allele %s is correct' % allele)
            return res

        return res

    @abstractmethod
    def read_result(self, temp):
        return None

    @abstractmethod
    def add_features(self, patient, file_tag_in, file_tag_out, write_res=True):
        return None

    @staticmethod
    def split_length(peptides):

        pept_len = [len(p) for p in peptides]

        dp = {}
        for pl in np.unique(pept_len):
            if pl > 7:
                dp[pl] = peptides[pept_len == pl]

        return dp

    @staticmethod
    def split_length_allele(peptides, alleles):

        assert(len(peptides) == len(alleles))

        pept_len_allele = list(set([(a, p, len(p)) for a, p in zip(alleles, peptides)]))

        d = {}
        for (a, p, l) in pept_len_allele:
            if len(p) > 7:
                if (a, l) in d.keys():
                    d[(a, l)] = d[(a, l)] + [p]
                else:
                    d[(a, l)] = [p]
        return d
