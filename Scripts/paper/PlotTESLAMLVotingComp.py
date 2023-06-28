import argparse
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns
import os

from Utils.GlobalParameters import GlobalParameters

parser = argparse.ArgumentParser(description='Plot comparison between TESLA groups and ML-Voting')

parser.add_argument('-fn', '--file_name', type=str, help='Name of plot output file')
parser.add_argument('-ft', '--file_type', type=str, default="pdf", choices=GlobalParameters.plot_file_formats,
                    help='File type for plot (png, svg or pdf)')
parser.add_argument('-pt', '--plot_type', type=str, choices=['TTIF_plot', 'FR_plot', 'AUPRC_plot'],
                    help='Type of plot corresponding to file_name')
parser.add_argument('-las', '--label_size', type=float, default=20.0, help='Axis label size')
parser.add_argument('-tis', '--tick_size', type=float, default=15.0, help='Tick size')
parser.add_argument('-fiw', '--figure_width', type=float, default=6.0, help='Figure width in inches')
parser.add_argument('-fih', '--figure_height', type=float, default=10.0, help='Figure height in inches')

if __name__ == "__main__":

    args = parser.parse_args()

    for arg in vars(args):
        print(arg, getattr(args, arg))

    file = GlobalParameters.tesla_result_file
    with open(file, 'rb') as file:
        data = pd.read_excel(file, sheet_name='auprc-by-team-patient-all', header=0)
        ML_Voting_results = {
            'TEAM': ['ML-Voting']*5,
            'PATIENT_ID': [1, 2, 3, 12, 16],
            'NUMBER RANKED': [0]*5,
            'TTIF': [0.3, 0.3, 0.25, 0.125, 0],
            'FR': [0.778, 1, 0.818, 0.667, 0.667],
            'AUPRC': [0.32, 0.273, 0.32, 0.094, 0.038],
            'total.validated': ['NA']*5
        }
        data = pd.concat([data, pd.DataFrame(ML_Voting_results)], axis=0, ignore_index=True)

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    pal = sns.color_palette("pastel", len(data.TEAM.unique()))
    color_dict = {team: "r" if team == "ML-Voting" else pal[i] for i, team in enumerate(data.TEAM.unique())}

    grouped = data.groupby('TEAM')

    if args.plot_type == 'TTIF_plot':
        fig = plt.figure()
        fig.set_figheight(args.figure_height)
        fig.set_figwidth(args.figure_height)
        data['TTIF-med'] = grouped[['TTIF']].transform(lambda x: x.median())
        data_sorted = data.sort_values(by=['TTIF-med', 'TEAM'], ascending=False)
        g = sns.boxplot(x="TTIF", y="TEAM", data=data_sorted, palette=color_dict)
        sns.swarmplot(x="TTIF", y="TEAM", data=data_sorted, color=".25")
        plt.xlabel('TTIF', fontsize=args.label_size)
        plt.ylabel('Team', fontsize=args.label_size)
        plt.xticks(fontsize=args.tick_size)
        plt.yticks(fontsize=args.tick_size)
        pdf_file = os.path.join(GlobalParameters.plot_dir, "{0}.{1}".format(args.file_name, args.file_type))
        plt.savefig(pdf_file, bbox_inches='tight', transparent=True)
        plt.close()

    if args.plot_type == 'FR_plot':
        fig = plt.figure()
        fig.set_figheight(args.figure_height)
        fig.set_figwidth(args.figure_height)
        data['FR-med'] = grouped[['FR']].transform(lambda x: x.median())
        data_sorted = data.sort_values(by=['FR-med', 'TEAM'], ascending=False)
        g = sns.boxplot(x="FR", y="TEAM", data=data_sorted, palette=color_dict)
        sns.swarmplot(x="FR", y="TEAM", data=data_sorted, color=".25")
        plt.xlabel('FR', fontsize=args.label_size)
        plt.ylabel('Team', fontsize=args.label_size)
        plt.xticks(fontsize=args.tick_size)
        plt.yticks(fontsize=args.tick_size)
        pdf_file = os.path.join(GlobalParameters.plot_dir, "{0}.{1}".format(args.file_name, args.file_type))
        plt.savefig(pdf_file, bbox_inches='tight', transparent=True)
        plt.close()

    if args.plot_type == 'AUPRC_plot':
        fig = plt.figure()
        fig.set_figheight(args.figure_height)
        fig.set_figwidth(args.figure_height)
        data['AUPRC-med'] = grouped[['AUPRC']].transform(lambda x: x.median())
        data_sorted = data.sort_values(by=['AUPRC-med', 'TEAM'], ascending=False)
        g = sns.boxplot(x="AUPRC", y="TEAM", data=data_sorted, palette=color_dict)
        sns.swarmplot(x="AUPRC", y="TEAM", data=data_sorted, color=".25")
        plt.xlabel('AUPRC', fontsize=args.label_size)
        plt.ylabel('Team', fontsize=args.label_size)
        plt.xticks(fontsize=args.tick_size)
        plt.yticks(fontsize=args.tick_size)
        pdf_file = os.path.join(GlobalParameters.plot_dir, "{0}.{1}".format(args.file_name, args.file_type))
        plt.savefig(pdf_file, bbox_inches='tight', transparent=True)
        plt.close()
