#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import pickle
from pathlib import Path
from tqdm.auto import tqdm
from collections import defaultdict
import logging

class Model_prediction:
    def __init__(self):
        """
        Initialize Model_prediction by loading models and required parameters.
        """
        # Configure logging
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

        # Define parameter and model paths
        base_dir = Path(__file__).parent / 'parameter'
        models_dir = Path(__file__).parent / 'models'

        hla_list_path = base_dir / 'hla_list.pkl'
        v_gene_list_path = base_dir / 'v_gene_list.pkl'
        models_1_path = models_dir / 'models_1.pkl'
        models_2_path = models_dir / 'models_2.pkl'

        # Load parameters and models
        self.hla_list = self.load_pickle(hla_list_path)
        self.v_gene_list = self.load_pickle(v_gene_list_path)
        models_1 = self.load_pickle(models_1_path)
        models_2 = self.load_pickle(models_2_path)

        # Merge both models
        models_1.update(models_2)
        self.models = models_1

    @staticmethod
    def load_pickle(file_path):
        """
        Load a pickle file with error handling.

        :param file_path: Path to the pickle file
        :return: Loaded object
        """
        try:
            with open(file_path, 'rb') as file:
                return pickle.load(file)
        except FileNotFoundError:
            logging.error(f"File not found: {file_path}")
            raise
        except Exception as e:
            logging.error(f"Failed to load {file_path}: {e}")
            raise

    def validate_v_genes(self, sample_v_gene_list):
        """
        Validate V genes against the acceptable list.

        :param sample_v_gene_list: List of V genes from input data
        """
        invalid_v_genes = [v for v in sample_v_gene_list if v not in self.v_gene_list]
        if invalid_v_genes:
            error_message = f"Invalid V genes detected: {', '.join(invalid_v_genes)}"
            logging.error(error_message)
            raise ValueError(error_message)

    def Get_prediction(self, input_df):
        """
        Generate predictions for each sample in the input DataFrame.

        :param input_df: DataFrame containing 'sample', 'cdr3', and 'v_gene'
        :return: Dictionary with sample-wise prediction scores
        """
        # Validate V genes
        sample_v_gene_list = list(set(input_df['v_gene'].tolist()))
        self.validate_v_genes(sample_v_gene_list)

        # Create TCR feature by combining 'cdr3' and 'v_gene'
        input_df['tcr'] = input_df['cdr3'] + input_df['v_gene']
        tcr_dict = input_df.groupby('sample')['tcr'].apply(list).to_dict()

        # Initialize the result dictionary
        sample_hit_rank = {}

        # Process each sample with progress feedback
        for sample_name, tcr_list in tqdm(tcr_dict.items(), desc="Predicting Samples"):
            tmp_cdr3_set = set(tcr_list)
            sample_hit = {}

            # Iterate through each HLA model
            for hla in tqdm(sorted(self.hla_list), desc=f"Processing {sample_name}", leave=False):
                clf = self.models[hla]

                # Initialize feature matrix
                tmp_df = pd.DataFrame(0, index=[sample_name], columns=clf.feature_names_in_)
                intersect_features = tmp_cdr3_set.intersection(clf.feature_names_in_)

                # Set intersecting features to 1
                if intersect_features:
                    tmp_df.loc[sample_name, list(intersect_features)] = 1

                # Predict probability
                y_pred = clf.predict_proba(tmp_df)
                sample_hit[hla] = float(y_pred[:, 1])

            # Store the result for the sample
            sample_hit_rank[sample_name] = sample_hit

        return sample_hit_rank
