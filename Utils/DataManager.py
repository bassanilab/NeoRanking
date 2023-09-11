import random
import numpy as np
import pandas as pd
from pandas.core.dtypes.concat import union_categoricals

from DataWrangling.CatEncoder import CatEncoder
from DataWrangling.MLRowSelection import MLRowSelection
from DataWrangling.DataTransformer import DataTransformer
from Utils.GlobalParameters import *
from Utils.Util_fct import *


class DataManager:
    """
    Class to load neo-peptide or mutation data from tab files
    Attributes:
        _original_data_dict (dict): dictionary containing non-processed data
        _ml_selected_data_dict (dict): data with rows used for training classifiers (e.g. only SNV rows)
        _processed_data_dict (dict): dictionary containing processed (normalized, imputed) data
        _immunogenic_patients (dict): dictionary with patients that contain immunogenic peptides
        _allotypes (pd.DataFrame): dataframe with patient HLA-I allotype info
        _cat_encoders (dict): encoders for categorical values
    """

    _original_data_dict: dict = {}
    _ml_selected_data_dict: dict = {}
    _processed_data_dict: dict = {}
    _immunogenic_patients: dict = {'mutation': None, 'neopep': None}
    _allotypes: pd.DataFrame = \
        pd.read_csv(filepath_or_buffer=GlobalParameters.hlaI_allele_file, sep="\t", header=0,
                    dtype={'Patient': 'string', 'Alleles': 'string'})
    _cat_encoders: dict = None

    @staticmethod
    def _get_filtered_data_index(data: pd.DataFrame, patient: str, dataset: str, response_types: list) -> list:
        """
        Private method. Retrieve boolean array indicating which rows are selected
        Args:
            data (pd.DataFrame): input data matrix
            patient (str): rows with mutations/neo-peptides of this patient will be set to True (if patient is not None)
            dataset (str): rows with mutations/neo-peptides of this dataset will be set to True (if dataset is not None)
            response_types (list): rows with one of these response types will be set to True (if response_types not empty)
        Returns:
            boolean array indicating which rows are selected
        """

        idx = np.full(data.shape[0], True)
        if patient != "":
            idx = np.logical_and(idx, data['patient'] == patient)
        elif dataset != "":
            if dataset == 'NCI_train':
                idx = np.logical_and(idx, (data['dataset'] == 'NCI') & (data['train_test'] == 'train'))
            elif dataset == 'NCI_test':
                idx = np.logical_and(idx, (data['dataset'] == 'NCI') & (data['train_test'] == 'test'))
            else:
                idx = np.logical_and(idx, data['dataset'] == dataset)

        response_types = set(response_types)
        if 0 < len(response_types) < 3:
            idx = np.logical_and(idx, data.response_type.apply(lambda row: row in response_types))

        return idx

    @staticmethod
    def load_filter_data(peptide_type: str, patient: str = "", dataset: str = "",
                         response_types: list = GlobalParameters.response_types,
                         ml_row_selection: bool = True) -> pd.DataFrame:
        """
        Function that returns the data matrix for neo-peptides or mutations. The data matrix can be filtered by
        patient, dataset, and response_type. The original complete data matrix is kept in memory for faster
        future access.

        Args:
            peptide_type (str): either 'neopep' or 'mutation'
            patient (str, optional): patient id. if not provided all patients are considered
            dataset (str, optional): dataset id (NCI, NCI_train, NCI_test, TESLA, HiTIDE).
                                     if not provided all datasets are considered
            response_types (list, optional): response_types ['CD8', 'negative', 'not_tested'] included in the data matrix.
                                             if empty all response types are considered
            ml_row_selection (bool): if True only rows used for ML are retrieved, if False all rows are retrieved

        Returns:
            Returns the dataframe selected according to the function arguments.
        """
        assert not peptide_type or peptide_type in GlobalParameters.peptide_types, \
            "DataManager.get_original_data: Unknown peptide_type."
        assert not dataset or dataset in GlobalParameters.datasets, \
            "DataManager.get_original_data: Unknown peptide_type."
        assert len(response_types) > 0 and all([rt in GlobalParameters.response_types for rt in response_types]), \
            "DataManager.get_original_data: Unknown response_type."

        if ml_row_selection:
            data = DataManager.load_ml_selected_data(peptide_type=peptide_type)
        else:
            data = DataManager.load_original_data(peptide_type=peptide_type)

        idx = DataManager._get_filtered_data_index(data=data, patient=patient, dataset=dataset,
                                                   response_types=response_types)
        return data.loc[idx, :]

    @staticmethod
    def load_original_data(peptide_type: str) -> pd.DataFrame:
        """
        Function that returns the unprocessed, unfiltered data matrix for neo-peptides or mutations. The original
        complete data matrix is kept in memory for faster future access.

        Args:
            peptide_type (str): either 'neopep' or 'mutation'

        Returns:
            Returns the original, unprocessed neo-peptide or mutation dataframe .
        """
        if peptide_type not in DataManager._original_data_dict:
            # in case data is not already loaded
            if peptide_type == 'neopep':
                data_file = GlobalParameters.neopep_data_org_file
                assert os.path.isfile(data_file), "No peopep data file {0}.".\
                    format(data_file)
                data = pd.read_csv(data_file, sep="\t", header=0)
                data = data.astype(dtype=GlobalParameters.feature_types_neopep, errors='ignore')
                data.replace('', 'nan', inplace=True)
            elif peptide_type == 'mutation':
                data_file = GlobalParameters.mutation_data_org_file
                assert os.path.isfile(data_file), \
                    "DataManager.load_original_data: No mutation data file {0}.".\
                    format(data_file)
                data = pd.read_csv(data_file, sep="\t", header=0)
                data = data.astype(dtype=GlobalParameters.feature_types_mutation, errors='ignore')
                data.replace('', 'nan', inplace=True)

            DataManager._original_data_dict[peptide_type] = data
        else:
            data = DataManager._original_data_dict[peptide_type]

        return data

    @staticmethod
    def load_ml_selected_data(peptide_type: str) -> pd.DataFrame:
        """
        Function that returns the data matrix for neo-peptides or mutations only containing the rows used for ML.

        Args:
            peptide_type (str): either 'neopep' or 'mutation'

        Returns:
            Returns the data matrix for neo-peptides or mutations only containing the rows used for ML.
        """
        if peptide_type not in DataManager._ml_selected_data_dict:
            # in case data is not already loaded
            if peptide_type == 'neopep':
                data_file = GlobalParameters.neopep_data_ml_sel_file
                assert os.path.isfile(data_file), "No peopep data file {0}. Download it and configure GlobalParameters.py".\
                    format(data_file)
                data = pd.read_csv(data_file, sep="\t", header=0)
                data = data.astype(dtype=GlobalParameters.feature_types_neopep, errors='ignore')
                data.replace('', 'nan', inplace=True)
            elif peptide_type == 'mutation':
                data_file = GlobalParameters.mutation_data_ml_sel_file
                assert os.path.isfile(data_file), \
                    "DataManager.load_original_data: No mutation data file {0}. Download it and configure GlobalParameters.py".\
                    format(data_file)
                data = pd.read_csv(data_file, sep="\t", header=0)
                data = data.astype(dtype=GlobalParameters.feature_types_mutation, errors='ignore')
                data.replace('', 'nan', inplace=True)

            DataManager._ml_selected_data_dict[peptide_type] = data
        else:
            data = DataManager._ml_selected_data_dict[peptide_type]

        return data

    @staticmethod
    def load_processed_data(peptide_type: str, objective: str) -> pd.DataFrame:
        """
        Function that returns the data matrix for neo-peptides or mutations only containing the rows used for ML, and
        processed either for plotting or classifier training/testing. For plotting a normalization transformation can
        be set for each feature, whereas for ML one global normalization transformation is set.

        Args:
            peptide_type (str): either 'neopep' or 'mutation'
            objective (str): 'ml' for machine learning, 'plot' for plotting

        Returns:
            Returns the data matrix for neo-peptides or mutations only containing the rows used for ML.
        """
        if peptide_type not in DataManager._processed_data_dict or \
                objective not in DataManager._processed_data_dict[peptide_type]:
            # in case data is not already loaded
            ml_sel_data_file_name, norm_data_file_name = DataManager._get_processed_data_files(peptide_type, objective)
            assert os.path.isfile(norm_data_file_name), "No data file {}. Use NormalizeData.py to create one.".\
                format(norm_data_file_name)
            X_ = pd.read_csv(norm_data_file_name, sep="\t", header=0, dtype=get_processed_types(peptide_type, objective))
            if peptide_type == 'neopep':
                X_ = X_.loc[:, GlobalParameters.ml_features_neopep]
                X_.index = X_.index.astype('int')
            elif peptide_type == 'mutation':
                X_ = X_.loc[:, GlobalParameters.ml_features_mutation]
                X_.index = X_.index.astype('int')

            DataManager._processed_data_dict[peptide_type] = {}
            DataManager._processed_data_dict[peptide_type][objective] = X_
        else:
            X_ = DataManager._processed_data_dict[peptide_type][objective]

        return X_

    @staticmethod
    def filter_processed_data(peptide_type: str, objective: str, patient: str = "", dataset: str = "",
                              response_types: list = GlobalParameters.response_types, sample: bool = True) -> tuple:
        """
        Function that returns the data matrix for neo-peptides or mutations after normalization and missing value
        imputation. The data matrix can be filtered by patient, dataset, and response_type. The original complete
        data matrix is kept in memory for faster future access.

        Args:
            peptide_type (str): either 'neopep' or 'mutation'
            objective (str): either machine learning ('ml') or plotting ('plot')
            patient (str, optional): patient id. if not provided all patients are considered
            dataset (bool, optional): dataset id. if not provided all datasets are considered
            response_types (list, optional): response_types ['CD8', 'negative', 'not_tested'] to be included in the
                                             data matrix.
            sample (bool): if true N rows with response_type != 'CD8' are randomly sampled

        Returns:
            Returns the data row-filtered original matrix, row-filtered normalized matrix, and vector indicating
            immunogenicity
        """
        peptide_type = peptide_type.lower()
        data = DataManager.load_ml_selected_data(peptide_type=peptide_type)
        X = DataManager.load_processed_data(peptide_type=peptide_type, objective=objective)
        y = np.array(data.response_type.apply(lambda rt: int(rt == 'CD8')), dtype=int)

        idx = DataManager._get_filtered_data_index(data=data, patient=patient, dataset=dataset,
                                                   response_types=response_types)

        if not all(idx):
            data = data.loc[idx, :]
            X = X.loc[idx, :]
            y = y[idx]

        if sample:
            data, X, y = DataManager._sample_rows(data=data, X=X, y=y)

        return data, X, y

    @staticmethod
    def filter_selected_data(peptide_type: str, patient: str = "", dataset: str = "",
                             response_types: list = GlobalParameters.response_types) -> pd.DataFrame:
        """
        Returns the data matrix for neo-peptides or mutations after normalization and missing value
        imputation. The data matrix can be filtered by patient, dataset, and response_type. The original complete
        data matrix is kept in memory for faster future access.

        Args:
            peptide_type (str): either 'neopep' or 'mutation'
            patient (str, optional): patient id. if not provided all patients are considered
            dataset (str, optional): dataset id. if not provided all datasets are considered
            response_types (list, optional): response_types ['CD8', 'negative', 'not_tested'] included in the data matrix.

        Returns:
            Returns the data row-filtered matrix.
        """
        peptide_type = peptide_type.lower()
        data = DataManager.load_ml_selected_data(peptide_type=peptide_type)
        idx = DataManager._get_filtered_data_index(data=data, patient=patient, dataset=dataset,
                                                   response_types=response_types)

        if not all(idx):
            data = data.loc[idx, :]

        return data

    @staticmethod
    def _combine_categories(df1, df2) -> tuple:
        for c in df1.columns:
            if df1[c].dtype.name == 'category':
                uc = union_categoricals([df1[c], df2[c]])
                df1.loc[:, c] = pd.Categorical(df1[c], categories=uc.categories)
                df2.loc[:, c] = pd.Categorical(df2[c], categories=uc.categories)

        return df1, df2

    @staticmethod
    def _sample_rows(data, X, y) -> tuple:
        if sum(y == 0) < GlobalParameters.nr_non_immuno_neopeps:
            return data, X, y

        idx = random.sample(range(sum(y == 0)), GlobalParameters.nr_non_immuno_neopeps)
        X_1 = X.loc[y == 1, :]
        X_0 = X.loc[y == 0, :]
        if X_0.shape[0] > GlobalParameters.nr_non_immuno_neopeps:
            X_0 = X_0.iloc[idx, :]
        X_s = pd.concat([X_1, X_0])

        data_1 = data.loc[y == 1, :]
        data_0 = data.loc[y == 0, :]
        if data_0.shape[0] > GlobalParameters.nr_non_immuno_neopeps:
            data_0 = data_0.iloc[idx, :]
        data_s = pd.concat([data_1, data_0])

        y_0 = y[y == 0]
        y_1 = y[y == 1]
        y_0 = y_0[idx]
        y_s = np.append(y_1, y_0)

        return data_s, X_s, y_s

    @staticmethod
    def has_immunogenic_peptides(peptide_type: str, patient: str) -> bool:
        """
        Checks whether patient has immunogenic mutations or neo-peptides

        Args:
            peptide_type (str): either 'neopep' or 'mutation'
            patient (str, optional): patient id.

        Returns:
            True if patient has immunogenic mutations or neo-peptides, False otherwise
        """
        peptide_type = peptide_type.lower()
        if not DataManager._immunogenic_patients[peptide_type]:
            data = DataManager.load_original_data(peptide_type=peptide_type)
            patients = data['patient'].unique()
            DataManager._immunogenic_patients[peptide_type] = []
            for p in patients:
                data_p = DataManager.load_filter_data(peptide_type=peptide_type, patient=p)
                if sum(data_p['response_type'] == 'CD8') > 0:
                    DataManager._immunogenic_patients[peptide_type].append(p)

        return patient in DataManager._immunogenic_patients[peptide_type]

    @staticmethod
    def _transform_data_(peptide_type: str, data_transformer: DataTransformer) \
            -> pd.DataFrame:
        data = DataManager.load_ml_selected_data(peptide_type=peptide_type)
        patients = data['patient'].unique()
        DataManager._immunogenic_patients[peptide_type] = []
        for i, p in enumerate(patients):
            print("processing patient {0}".format(p))
            data_p = \
                DataManager.load_filter_data(peptide_type=peptide_type, patient=p, ml_row_selection=True)
            data_p, X_p, y_p = data_transformer.apply(data_p)
            if i == 0:
                combined_df = data_p
                combined_X = X_p
                combined_y = y_p
            else:
                combined_df, data_p = DataManager._combine_categories(combined_df, data_p)
                combined_df = pd.concat([combined_df, data_p], ignore_index=True)
                if X_p is not None:
                    combined_X, X = DataManager._combine_categories(combined_X, X_p)
                    combined_X = pd.concat([combined_X, X_p], ignore_index=True)
                combined_y = np.append(combined_y, y_p)

        return combined_df, combined_X, combined_y

    @staticmethod
    def transform_data(peptide_type: str, dataset: str, objective: str) -> None:
        """
        Transforms data by encoding categorical values, normalizing numerical values, and imputing missing values.
        These operations are performed separately for each patient in all datasets. Stores the transformed data matrix
        in tab separated text files.

        Args:
            peptide_type (str): either 'neopep' or 'mutation'
            dataset (bool, optional): dataset id used to train the encoding of categorical features
            objective (str): either machine learning ('ml') or plotting ('plot')

        """
        data_transformer = DataTransformer(peptide_type, objective, dataset, DataTransformer.get_normalizer(objective))
        data, X, y = DataManager._transform_data_(peptide_type=peptide_type, data_transformer=data_transformer)
        ml_sel_data_file_name, norm_data_file_name = DataManager._get_processed_data_files(peptide_type, objective)
        X.to_csv(norm_data_file_name, sep='\t', header=True, index=False)
        data.to_csv(ml_sel_data_file_name, sep='\t', header=True, index=False)

    @staticmethod
    def select_ml_data(peptide_type: str) -> None:
        """
        Selects the rows used for ML. Stores the filtered data matrix in tab separated text files.

        Args:
            peptide_type (str): either 'neopep' or 'mutation'

        """
        data = DataManager.load_original_data(peptide_type=peptide_type)
        data = MLRowSelection.apply(data=data, peptide_type=peptide_type)
        ml_sel_data_file_name = DataManager._get_processed_data_files(peptide_type, 'sel')[0]
        data.to_csv(ml_sel_data_file_name, sep='\t', header=True, index=False)

    @staticmethod
    def _get_processed_data_files(peptide_type: str, objective: str = 'sel') -> tuple:
        if peptide_type == 'neopep' and objective == 'ml':
            return GlobalParameters.neopep_data_ml_sel_file, GlobalParameters.neopep_data_ml_file
        elif peptide_type == 'neopep' and objective == 'plot':
            return GlobalParameters.neopep_data_ml_sel_file, GlobalParameters.neopep_data_plot_file
        elif peptide_type == 'neopep' and objective == 'sel':
            return GlobalParameters.neopep_data_ml_sel_file, None
        elif peptide_type == 'mutation' and objective == 'ml':
            return GlobalParameters.mutation_data_ml_sel_file, GlobalParameters.mutation_data_ml_file
        elif peptide_type == 'mutation' and objective == 'plot':
            return GlobalParameters.mutation_data_ml_sel_file, GlobalParameters.mutation_data_plot_file
        elif peptide_type == 'mutation' and objective == 'sel':
            return GlobalParameters.mutation_data_ml_sel_file, None
        else:
            return None, None

    @staticmethod
    def get_classI_allotypes(patient: str) -> list:
        """
        Retrieves list of HLA-I allotypes for a patient
        Args:
            patient (str): patient id

        Returns:
            list of HLA-I allotypes
        """

        if any(DataManager._allotypes['Patient'].str.contains(patient)):
            a = DataManager._allotypes.loc[DataManager._allotypes['Patient'] == patient, 'Alleles'].iloc[0]
            return a.split(sep=",")
        else:
            return []

    @staticmethod
    def get_categorical_feature_idx(peptide_type: str, x_: pd.DataFrame) -> list:
        """
        Retrieves indexes of columns with categorical features
        Args:
            peptide_type (str):  either 'neopep' or 'mutation'
            x_ (pd.DataFrame): ML dataframe with feature columns in order
        Returns:
            list of indexes of columns with categorical features
        """
        if peptide_type == 'neopep':
            idx = [i for i, c in enumerate(x_.columns) if GlobalParameters.feature_types_neopep[c] == 'category']
        else:
            idx = [i for i, c in enumerate(x_.columns) if GlobalParameters.feature_types_mutation[c] == 'category']

        return idx

    @staticmethod
    def get_category_cnts(dataset: str, peptide_type: str, x_: pd.DataFrame) -> dict:
        """
        Counts number of categories in categorical features in x_
        Args:
            dataset: dataset id used to train the encoding of categorical features
            peptide_type (str):  either 'neopep' or 'mutation'
            x_ (pd.DataFrame): ML dataframe with feature columns in order
        Returns:
            dictionary of with features and counts
        """
        if DataManager._cat_encoders is None:
            DataManager._cat_encoders = CatEncoder.read_cat_encodings(dataset=dataset, peptide_type=peptide_type)

        if peptide_type == 'neopep':
            cat_cnts = {c: DataManager._cat_encoders[c].get_nr_classes() for c in x_.columns
                        if GlobalParameters.feature_types_neopep[c] == 'category'}
        else:
            cat_cnts = {c: DataManager._cat_encoders[c].get_nr_classes() for c in x_.columns
                        if GlobalParameters.feature_types_mutation[c] == 'category'}

        return cat_cnts

    @staticmethod
    def _create_mutation_id(data_row):
        return data_row['chromosome'] + "/" + data_row['genomic_coord'] + "/" + data_row['ref'] + "/" + data_row['alt']
