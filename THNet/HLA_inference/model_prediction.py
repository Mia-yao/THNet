#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import pickle
import sys, os
from tqdm.auto import tqdm
from collections import defaultdict

class Model_prediction:
    def __init__(self):
        """
        Initialize Model_prediction.

        :param models: models for each HLA inference model
        """
        hla_list_path = os.path.join(os.path.dirname(__file__), 'parameter', 'hla_list.pkl')
        Vgene_list_path = os.path.join(os.path.dirname(__file__), 'parameter', 'v_gene_list.pkl')
        models_1_path = os.path.join(os.path.dirname(__file__), 'models', 'models_1.pkl')
        models_2_path = os.path.join(os.path.dirname(__file__), 'models', 'models_2.pkl')

        with open(hla_list_path, 'rb') as file:
            hla_list = pickle.load(file)
        with open(Vgene_list_path, 'rb') as file:
            v_gene_list = pickle.load(file)
        with open(models_1_path, 'rb') as file:
            models = pickle.load(file)      
        with open(models_2_path, 'rb') as file:
            models_2 = pickle.load(file) 
        models.update(models_2)
        
        self.models=models
        self.hla_list=hla_list
        self.v_gene_list=v_gene_list
        
    def Get_prediction(self,input_df):
        # Convert the dataframe into a dictionary
        tcr_dict = defaultdict(list)
        sample_v_gene_list=list(set(input_df['v_gene'].tolist()))
        
        failed_v=[]
        for v_gene in sample_v_gene_list:
            if v_gene not in self.v_gene_list:
                failed_v.append(v_gene)
        
        if len(failed_v) != 0:
            print(f"{','.join(failed_v)} is/are not in the acceptable V gene list.")
            sys.exit()  # Exit the program if there are invalid V genes
        
        # If all V genes are valid, proceed with the prediction
        input_df['tcr'] = input_df['cdr3'] + input_df['v_gene']
        tcr_dict = input_df.groupby('sample')['tcr'].apply(list).to_dict()
        
        # Initialize the result dictionary
        sample_hit_rank = {}

        # Process samples in a single thread
        for sample_name in tqdm(tcr_dict.keys()):
            tmp_cdr3_set = set(tcr_dict[sample_name])
            sample_hit = {}

            for hla in sorted(self.hla_list):
                # Get the model and non-zero features
                clf = self.models[hla]

                # Initialize the feature matrix
                tmp_df = pd.DataFrame(0, index=[sample_name], columns=clf.feature_names_in_)

                # Set the value of the intersecting features to 1
                tmp_df.loc[sample_name, list(tmp_cdr3_set.intersection(clf.feature_names_in_))] = 1

                # Model prediction
                y_pred = clf.predict_proba(tmp_df)
                sample_hit[hla] = float(y_pred[:, 1])

            # Save the prediction results for the current sample
            sample_hit_rank[sample_name] = sample_hit
        
        return sample_hit_rank
