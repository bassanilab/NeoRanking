import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import os

from Utils.Parameters import Parameters

file = Parameters().get_tesla_result_file()
with open(file, 'rb') as file:
    data = pd.read_excel(file, sheet_name='auprc-by-team-patient-all', header=0)

grouped = data.groupby('TEAM')

fig = plt.figure()
fig.set_figheight(10)
fig.set_figwidth(6)
data['TTIF-med'] = grouped[['TTIF']].transform(lambda x: x.median())
data_sorted = data.sort_values(by=['TTIF-med', 'TEAM'], ascending=False)
g = sns.boxplot(x="TTIF", y="TEAM", data=data_sorted)
sns.swarmplot(x="TTIF", y="TEAM", data=data_sorted, color=".25")
plt.xlabel('TTIF', fontsize=20)
plt.ylabel('Team', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
pdf_file = os.path.join(Parameters().get_plot_dir(), "TESLA_comp_TTIF.pdf")
plt.savefig(pdf_file, bbox_inches='tight')
plt.close()

fig = plt.figure()
fig.set_figheight(10)
fig.set_figwidth(6)
data['FR-med'] = grouped[['FR']].transform(lambda x: x.median())
data_sorted = data.sort_values(by=['FR-med', 'TEAM'], ascending=False)
g = sns.boxplot(x="FR", y="TEAM", data=data_sorted)
sns.swarmplot(x="FR", y="TEAM", data=data_sorted, color=".25")
plt.xlabel('FR', fontsize=20)
plt.ylabel('Team', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
pdf_file = os.path.join(Parameters().get_plot_dir(), "TESLA_comp_FR.pdf")
plt.savefig(pdf_file, bbox_inches='tight')
plt.close()

fig = plt.figure()
fig.set_figheight(10)
fig.set_figwidth(6)
data['AUPRC-med'] = grouped[['AUPRC']].transform(lambda x: x.median())
data_sorted = data.sort_values(by=['AUPRC-med', 'TEAM'], ascending=False)
g = sns.boxplot(x="AUPRC", y="TEAM", data=data_sorted)
sns.swarmplot(x="AUPRC", y="TEAM", data=data_sorted, color=".25")
plt.xlabel('AUPRC', fontsize=20)
plt.ylabel('Team', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
pdf_file = os.path.join(Parameters().get_plot_dir(), "TESLA_comp_AUPRC.pdf")
plt.savefig(pdf_file, bbox_inches='tight')
plt.close()
