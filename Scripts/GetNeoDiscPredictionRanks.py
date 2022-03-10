import argparse
from DataWrangling.DataLoader import *
from Utils.Util_fct import *


parser = argparse.ArgumentParser(description='Plot and test difference between immunogenic and non immunogenic feature'
                                             'values')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-m', '--method', type=str, default='selection',
                    help='Method for NeoDisc prioritization  (selection or prioritization)')


args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

data_loader = DataLoader(transformer=None, normalizer=None, features=None,
                         mutation_types=args.mutation_types, response_types=args.response_types,
                         immunogenic=['CD8', 'CD4/CD8'], min_nr_immono=0)


patients = get_valid_patients("TESLA/HiTIDE")
for p in patients:
    data, X, y = data_loader.load_patients(p, "rt", args.peptide_type)
    if data is None or data.shape[0] == 0:
        print(f"Patient p: {[]}")
        continue

    if args.method == 'selection':
        data_file = glob.glob(os.path.join(Parameters().get_data_dir(), 'neodisc_selection', p + '-*_Selection_CI.txt'))
        if len(data_file) == 1:
            data_file = data_file[0]
            neodisc_ordered = pd.read_csv(data_file, header=0, sep="\t")
        else:
            print(f"No or many files found for patient {p}")
            continue
    else:
        neodisc_ordered = data

    CD8_peptides = data.loc[data.apply(lambda row: row['response_type'] in ['CD8', 'CD4/CD8'], axis=1), 'mutant_seq']
    ranks = []
    for pept in CD8_peptides:
        idx = np.where(neodisc_ordered['mutant_seq'] == pept)
        if len(idx[0]) > 0:
            ranks.append(idx[0][0])
        else:
            print(f"Peptide {pept} not found in patient {p}")

    print(f"Patient {p}: {np.sort(ranks)}")





