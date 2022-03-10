import tempfile
from Utils import *
import subprocess
from Neoantigen import Neoantigen
from Aligner import Aligner
from Neoantigen import Neoantigen

def run_blastp(blaster, netMHC_data, out_file):

    ids = np.array(netMHC_data['Identity'])
    cores = np.array(netMHC_data['Icore'])
    pep_file = tempfile.mktemp() + '.fa'
    with open(pep_file, 'w') as f:
        for i in np.arange(netMHC_data.shape[0]):
            f.write('>' + ids[i] + '\n' + cores[i] + '\n')
    f.close()

    blaster.run(pep_file, out_file=out_file)

    os.remove(pep_file)

    return


def add_foreigness_info(netmhc_res, mut_peptides, wt_peptides, mut_blastp_outfile, wt_blastp_outfile):

    aligner = Aligner()
    aligner.readAllBlastAlignments(mut_blastp_outfile)
    aligner.computeR()

    for i in netmhc_res.index.values:
        neoantigen = Neoantigen([netmhc_res['wt']])

    return


class BlastRunner:

    def __init__(self, path=None, target_file=None, e_value=10, gap_open=11, gap_extend=1):

        self.blastpPath = find_exe(path=path,exe='blastp')
        self.target_file = target_file
        self.e_value = e_value
        self.gap_open = gap_open
        self.gap_extend = gap_extend

        if not os.seq_logo_dir.isfile(target_file):
            print('File %s does not exist.' %target_file)
            self.target_file = None
        return

    def run(self, query_file, out_file, show_cmd=True):

        if self.target_file is None:
            print('No valid target file for blastp. Abort')
            self.target_file = None
            return

        'blastp -query cow.small.faa -db human.1.protein.faa -out cow_vs_human_blast_results.txt'
        cmd = '%s -query %s -db %s -evalue %f -gapopen %d -gapextend %d -outfmt 5 -out %s' % \
              (self.blastpPath, query_file, self.target_file, self.e_value, self.gap_open, self.gap_extend, out_file)

        if show_cmd is True:
            print (cmd)

        try:
            temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
        except Exception as e:
            print(e)
            return


def main():
    return


if __name__ == "__main__":
    main()
