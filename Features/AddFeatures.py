import argparse

from Features.NetMHC.NetMHCpanIRunner import *
from Features.NetMHC.NetMHCstabIRunner import *
from Features.NetMHC.NetChopRunner import *
from Features.NetMHC.NetTAPRunner import *
from Features.BindingInfo.MutationAAChange import *
from Features.BindingInfo.CalcBindingAffinityScoreDists import *
from DataWrangling.RosenbergImmunogenicityAnnotatorLong import *

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-p', '--patient', type=str, nargs='+', help='prefix for patient id')
parser.add_argument('-t', '--input_file_tag', type=str, nargs='?', default='',
                    help='File tage for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('--netmhcpan', action='store_true', help='Run NetMHCpan 4.1')
parser.add_argument('--netmhcstab', action='store_true', help='Run NetMHCstabpan 1.0')
parser.add_argument('--netchop', action='store_true', help='Run NetChop 3.1')
parser.add_argument('--nettap', action='store_true', help='Run NetChop 3.1')
parser.add_argument('--seqlogo', action='store_true', help='Run SequenceLogo prediction')
parser.add_argument('--tcr', action='store_true', help='Run MutationAAChange prediction')
parser.add_argument('--prop', action='store_true', help='Run allele propensity prediction')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

mgr = DataManager()
parameters = Parameters()

patients = mgr.get_valid_patients()
patients = [p for p in patients if any([ap == p for ap in args.patient])]

print("Running patients {0}".format(patients))
for patient in patients:

    if patient in mgr.get_valid_patients():
        print("add features to: "+str(patient))
        tag = str(args.input_file_tag)

        if args.netmhcpan and args.peptide_type == 'long':
            netMHCpanRunner = NetMHCpanIRunner(exe_base_dir=parameters.get_exe_dir(), max_rank_peptide=5)
            print("Run NetMHCpan for patient {0} with file tag {1}".format(patient, tag))
            tag_out = tag+"_netmhc" if tag else "netmhc"
            netMHCpanRunner.add_features(patient, tag, tag_out, write_res=True)
            tag = tag_out

        if args.netmhcstab:
            netMHCstabRunner = NetMHCstabIRunner(exe_base_dir=parameters.get_exe_dir())
            print("Run NetMHCstabpan for patient {0} with file tag {1}".format(patient, tag))
            tag_out = tag+"_stab" if tag else "stab"
            netMHCstabRunner.add_features(patient, tag, tag_out, args.peptide_type, write_res=True)
            tag = tag_out

        if args.netchop:
            netchopRunner = NetChopRunner(exe_path=parameters.get_exe_dir())
            print("Run NetChop for patient {0} with file tag {1}".format(patient, tag))
            tag_out = tag+"_chop" if tag else "chop"
            netchopRunner.add_features(patient, tag, tag_out, args.peptide_type, write_res=True)
            tag = tag_out

        if args.nettap:
            nettapRunner = NetTAPRunner(exe_path=os.path.join(parameters.get_exe_dir(), 'netCTLpan-1.1', 'Linux_x86_64'))
            print("Run NetChop for patient {0} with file tag {1}".format(patient, tag))
            tag_out = tag+"_tap" if tag else "tap"
            nettapRunner.add_features(patient, tag, tag_out, args.peptide_type, write_res=True)
            tag = tag_out

        if args.seqlogo:
            sequenceLogoMgr = SequenceLogoMgr(seq_logo_dir=parameters.get_data_dir(), binding_threshold=20)
            print("Run SequenceLogo for patient {0} with file tag {1}".format(patient, tag))
            tag_out = tag+"_mbp" if tag else "mbp"
            sequenceLogoMgr.add_features(patient, tag, tag_out, args.peptide_type, write_res=True)
            tag = tag_out

        if args.tcr:
            mutationAAChange = MutationAAChange()
            print("Run MutationAAChange for patient {0} with file tag {1}".format(patient, tag))
            tag_out = tag+"_tcr" if tag else "tcr"
            mutationAAChange.add_features(patient, tag, tag_out, args.peptide_type, write_res=True)
            tag = tag_out

        if args.prop:
            annotator = RosenbergImmunogenicityAnnotatorLong(mgr)
            calcAllelePropensity = CalcAllelePropensity(parameters.get_human_allele_score_dist_file(),
                                                        annotator.get_patients('parkhurst'), 'SNV', 'CD8,CD4/CD8',
                                                        'CD8,CD4/CD8,negative', 'rt_netmhc', 'long')

            print("Run CalcAllelePropensity for patient {0} with file tag {1}".format(patient, tag))
            tag_out = tag+"_prop" if tag else "prop"
            calcAllelePropensity.add_features(patient, tag, tag_out, args.peptide_type, write_res=True)
            tag = tag_out

    else:
        print(f"Patient {patient} does not have a data file. Skip.")
