import numpy as np
import pandas as pd

from collections import Counter

from sklearn.preprocessing import QuantileTransformer, StandardScaler, PowerTransformer, MinMaxScaler, \
    FunctionTransformer

from Utils.GlobalParameters import GlobalParameters
from Utils.Util_fct import *
from DataWrangling.CatEncoder import CatEncoder


class MLRowSelection:

    @staticmethod
    def apply(data: pd.DataFrame, peptide_type: str) -> pd.DataFrame:
        if data is None or data.shape[0] == 0:
            return None, None, None

        if peptide_type == 'mutation':
            return MLRowSelection.filter_rows_mutation(data)

        elif peptide_type == 'neopep':
            return MLRowSelection.filter_rows_neopep(data)
        else:
            return None

    @staticmethod
    def filter_rows_mutation(df):
        if df.shape[0] > 0:
            df = df.loc[df.mutation_type.apply(lambda r: r == 'SNV')]

        if df.shape[0] > 0 and GlobalParameters.max_netmhc_rank > 0:
            df = df.loc[df.mut_Rank_EL_0.apply(lambda r: r < GlobalParameters.max_netmhc_rank)]

        if df.shape[0] > 0 and len(GlobalParameters.excluded_genes) > 0:
            df = df.loc[df.gene.apply(lambda g: g not in GlobalParameters.excluded_genes)]

        df.reset_index(inplace=True, drop=True)

        return df

    @staticmethod
    def filter_rows_neopep(df):
        if df.shape[0] > 0:
            df = df.loc[df.mutation_type.apply(lambda r: r == 'SNV')]

        if df.shape[0] > 0 and GlobalParameters.max_netmhc_rank > 0:
            df = df.loc[df.mutant_rank_netMHCpan.apply(lambda r: r < GlobalParameters.max_netmhc_rank)]

        if df.shape[0] > 0 and len(GlobalParameters.excluded_genes) > 0:
            df = df.loc[df.gene.apply(lambda g: g not in GlobalParameters.excluded_genes)]

        df.reset_index(inplace=True, drop=True)

        return df

