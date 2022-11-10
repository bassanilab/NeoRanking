import argparse

from DataWrangling.DataLoader import *
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Print dataset statistics')

parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')

args = parser.parse_args()

output_file = os.path.join(Parameters().get_plot_dir(), "Dataset_patients_{0}.txt".format(args.peptide_type))
with open(output_file, 'w') as file:
    file.write("Patient group\tPatient count\tPatients\n")
    for patient_group in ['TESLA', 'HiTIDE', 'NCI', 'NCI_test', 'NCI_train']:
        patients = get_valid_patients(dataset=patient_group, peptide_type=args.peptide_type)
#        patients = get_valid_annotated_patients(dataset=patient_group, peptide_type=args.peptide_type)
        file.write("{0}\t{1}\t{2}\n".format(patient_group, len(patients), ",".join(sorted(patients))))
