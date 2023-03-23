from Utils.DataManager import *


class ReplaceIpMSDBScores:

    def __init__(self, file_tag_input, file_tag_output, peptide_type):
        """ """
        mgr = DataManager()
        patients = mgr.get_valid_patients(peptide_type)
        for p in patients:
            self.remove_ipmsdb_features(p, file_tag_input, file_tag_output, peptide_type, mgr)

        return

    @staticmethod
    def remove_ipmsdb_features(patient, file_tag_input, file_tag_output, peptide_type='long', mgr=DataManager()):
        parameters = Parameters()

        data = mgr.get_processed_data(patient, file_tag_input, peptide_type)
        if data is None:
            return

        data = data.loc[:, ~data.columns.str.startswith('best')]

        out_file = os.path.join(parameters.get_result_dir(), patient+"_"+peptide_type+"_"+file_tag_output+".txt")
        data.to_csv(path_or_buf=out_file, sep="\t", index=False, header=True)


if __name__ == "__main__":
    ReplaceIpMSDBScores(file_tag_input="rt_netmhc_stab_chop_tap_mbp_tcr",
                        file_tag_output="rt_netmhc_stab_chop_tap_mbp_tcr_ipmsdb", peptide_type='long')
    ReplaceIpMSDBScores(file_tag_input="rt_stab_chop_tap_mbp_tcr", file_tag_output="rt_stab_chop_tap_mbp_tcr_ipmsdb",
                        peptide_type='short')
