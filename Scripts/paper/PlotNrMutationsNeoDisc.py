import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from Utils.Parameters import *
from collections import Counter

pdf_file = os.path.join(Parameters().get_plot_dir(),  'GartnerNeodiscComp.pdf')
data_file = os.path.join(Parameters().get_plot_dir(),  'Compare_NeoDisc_Gartner.txt')

data = pd.read_csv(filepath_or_buffer=data_file, sep="\t", header=0)

with PdfPages(pdf_file) as pp:

    fig = plt.figure(figsize=(10, 10))
    fig.clf()
    g = sns.scatterplot(data=data, x="Gartner et al. mutations", y="NeoDisc mutations", alpha=0.7, s=100)
    plt.xlabel("Gartner et al. mutations", size=15)
    plt.ylabel("NeoDisc mutations", size=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    g.set_xscale('log')
    g.set_yscale('log')
    plt.plot([0, 1.05*data["Gartner et al. mutations"].max()],
             [0, 1.05*data["Gartner et al. mutations"].max()], color='r', alpha=0.7,)
    g.figure.tight_layout()
    pp.savefig()
    plt.close()

    fig = plt.figure(figsize=(10, 10))
    fig.clf()
    diff = data.loc[:, "Gartner et al. CD8+ immunogenic mutations screened with minigenes"] - \
           data.loc[:, "NeoDisc CD8+ immunogenic mutations screened with minigenes"]

    counter = Counter(diff)

    g = sns.barplot(x=list(counter.keys()), y=list(counter.values()))
    plt.xlabel("Difference CD8+ mutations screened with minigene", size=15)
    plt.ylabel("Counts", size=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    for k, v in counter.items():
        g.text(k, v+0.5, v, color='black', ha="center", fontsize=12)
    g.figure.tight_layout()
    pp.savefig()
    plt.close()
