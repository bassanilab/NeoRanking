from Features.BindingInfo.CalcBindingAffinityScoreDists import *
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns


mgr = DataManager()
proteins = mgr.get_ipMSDB_proteins(52)
fasta_file = Parameters().get_human_fasta_file()
#fasta_file = Parameters().get_virus_fasta_file()
#proteins = set()
promiscuity = CalcBindingAffinityScoreDists(seq_logo_dir=Parameters().get_data_dir(), seq_fasta_file=fasta_file,
                                            proteins=proteins)

alleles = sorted(mgr.get_all_classI_allotypes())
output_file = path.join(Parameters().get_plot_dir(), "allele_score_data_human.txt")
pdf_file = path.join(Parameters().get_plot_dir(), "allele_score_dists_human.pdf")

with open(output_file, 'w') as out_file:
    out_file.write(promiscuity.get_header()+"\n")
    with PdfPages(pdf_file) as pp:
        for a in alleles:
            print("PlotMixMHCScoreDists.py: processing {0} ...".format(a))
            scores = []
            lengths = []
            for peptide_length in [9,  10, 11, 12]:
                s = promiscuity.match_fasta(a, peptide_length)
                scores = np.append(scores, s)
                lengths = np.append(lengths, np.full(len(s), peptide_length))
                values = promiscuity.calc_dist_stats(alleles, s, peptide_length)
                out_file.write("\t".join([str(v) for v in values])+"\n")
                out_file.flush()

            df = pd.DataFrame({'scores': scores, 'length': lengths})
            bins = np.arange(start=-30, stop=15, step=0.2)
            g = sns.histplot(
                       data=df, x='scores', hue="length", bins=bins, legend=True, stat="density",
                       fill=True, common_norm=False, palette="viridis", alpha=.5, linewidth=0
                    )

            g.set_title("allele = {0}, length = {1}".format(a, peptide_length))
            g.figure.tight_layout()
            pp.savefig(g.figure)
            g.figure.clf()
