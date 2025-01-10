#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import pickle
import sys, os

class HLA_inference:
    def __init__(self):
        """
        Initialize HLAProcessor.

        :param hla_threshold_dict: Dictionary containing HLA threshold values
        :param hla_list: List of HLA names
        """
        # Load threshold file
        hla_auc_path = os.path.join(os.path.dirname(__file__), 'parameter', 'hla_auc.pkl')
        fscore_dict_path = os.path.join(os.path.dirname(__file__), 'parameter', 'fscore_dict.pkl')
        hla_threshold_path = os.path.join(os.path.dirname(__file__), 'parameter', 'hla_threshold_dict.pkl')
        
        with open(hla_auc_path, 'rb') as file:    
            hla_auc= pickle.load(file)
        with open(fscore_dict_path, 'rb') as file:    
            fscore_dict= pickle.load(file)
        with open(hla_threshold_path, 'rb') as file:
            hla_threshold_dict = pickle.load(file)
        hla_list=list(hla_threshold_dict.keys())
        
        self.hla_threshold_dict = hla_threshold_dict
        self.hla_list = hla_list
        self.hla_auc = hla_auc
        self.fscore_dict = fscore_dict

    def find_pred(self, start, end, prediction, max_bar, second_bar):
        """
        Compute the prediction results for a specified range.

        :param start: Start index
        :param end: End index
        :param prediction: List of prediction values
        :param max_bar: Threshold for the maximum value
        :param second_bar: Threshold for the second highest value
        :return: Prediction results for the specified range
        """
        # Extract threshold values and predictions for the specified range
        threshold_tmp = np.array(list(self.hla_threshold_dict.values())[start:end])
        section = np.array(prediction[start:end]) - threshold_tmp

        # Initialize the result array
        result = np.zeros_like(section)

        # Get the indices and values of the top 3 values
        top_indices = np.argsort(section)[-3:]
        top_values = section[top_indices]

        # Extract the maximum, second, and third values
        max_value = top_values[2]
        second_value = top_values[1]
        third_value = top_values[0]

        # Get the indices of the maximum and second highest values
        max_value_index = top_indices[2]
        second_value_index = top_indices[1]

        # Set the result based on the conditions
        if second_value > second_bar:
            result[[max_value_index, second_value_index]] = 1
        elif max_value > max_bar:
            result[max_value_index] = 1

        return result

    def convert_pred(self, prediction, max_bar=0.1, second_bar=0.1):
        """
        Convert prediction values to results.

        :param prediction: List of prediction values
        :param max_bar: Threshold for the maximum value
        :param second_bar: Threshold for the second highest value
        :return: Converted prediction results
        """
        A_part = self.find_pred(0, 32, prediction, max_bar, second_bar)
        B_part = self.find_pred(32, 83, prediction, max_bar, second_bar)
        C_part = self.find_pred(83, 109, prediction, max_bar, second_bar)
        DPA1_part = self.find_pred(109, 117, prediction, max_bar, second_bar)
        DPB1_part = self.find_pred(117, 138, prediction, max_bar, second_bar)
        DQA1_part = self.find_pred(138, 152, prediction, max_bar, second_bar)
        DQB1_part = self.find_pred(152, 169, prediction, max_bar, second_bar)
        DRB1_part = self.find_pred(169, 208, prediction, max_bar, second_bar)

        return np.concatenate((
            A_part, B_part, C_part, DPA1_part, DPB1_part, DQA1_part, DQB1_part, DRB1_part
        )).astype(int)

    def get_pred_dict(self, sample_hit_rank):
        """
        Process all samples and generate prediction results.

        :param sample_hit_rank: Sample ranking data
        :return: Dictionary of prediction results for samples
        """
        pred_dict = {}

        for sample, value in sample_hit_rank.items():
            prob_value = np.array(list(value.values()))
            indexes = np.where(self.convert_pred(prob_value) == 1)[0]
            pred_dict[sample] = [self.hla_list[index] for index in indexes]
        return pred_dict
    
    def hla_inference_df(self, sample_hit_rank):
        columns = ['Sample', 'HLA-A_1', 'HLA-A_2', 'HLA-B_1', 'HLA-B_2', 'HLA-C_1', 'HLA-C_2',
                   'HLA-DRB1_1', 'HLA-DRB1_2', 'HLA-DQA1_1', 'HLA-DQA1_2', 'HLA-DQB1_1', 'HLA-DQB1_2', 
                   'HLA-DPA1_1', 'HLA-DPA1_2', 'HLA-DPB1_1', 'HLA-DPB1_2']

        # Initialize a list to store row data for each sample
        rows = []

        # Populate HLA data by separating Class I and Class II HLAs
        pred_dict = self.get_pred_dict(sample_hit_rank)
        for sample, hla_list in pred_dict.items():
            class_I = ['NA'] * 6  # Initial Class I HLA columns
            class_II = ['NA'] * 10  # Initial Class II HLA columns

            # Fill Class I and Class II columns
            for hla in hla_list:
                if '-A' in hla:
                    if class_I[0] == 'NA':
                        class_I[0] = hla
                    else:
                        class_I[1] = hla
                elif '-B' in hla:
                    if class_I[2] == 'NA':
                        class_I[2] = hla
                    else:
                        class_I[3] = hla
                elif '-C' in hla:
                    if class_I[4] == 'NA':
                        class_I[4] = hla
                    else:
                        class_I[5] = hla
                elif 'DRB1' in hla:
                    if class_II[0] == 'NA':
                        class_II[0] = hla
                    else:
                        class_II[1] = hla
                elif 'DQA1' in hla:
                    if class_II[2] == 'NA':
                        class_II[2] = hla
                    else:
                        class_II[3] = hla
                elif 'DQB1' in hla:
                    if class_II[4] == 'NA':
                        class_II[4] = hla
                    else:
                        class_II[5] = hla
                elif 'DPA1' in hla:
                    if class_II[6] == 'NA':
                        class_II[6] = hla
                    else:
                        class_II[7] = hla
                elif 'DPB1' in hla:
                    if class_II[8] == 'NA':
                        class_II[8] = hla
                    else:
                        class_II[9] = hla

            rows.append([sample] + class_I + class_II)

        return pd.DataFrame(rows, columns=columns)
    def create_top_hla_df(self, sample_hit_rank, top_n=3):
        
        """
        Preprocess the sample_hit_rank dictionary by grouping HLA types and keeping top 3 values for each group.

        :param sample_hit_rank: Original dictionary with sample HLA rankings
        :return: Processed sample_hit_rank dictionary
        """
        
        sample_hit_rank_2 = {}

        # Group HLA types by splitting keys at '*'
        for sample, original_dict in sample_hit_rank.items():
            new_dict = {}
            for key, value in original_dict.items():
                new_key = key.split('*')[0] + '*'
                if new_key not in new_dict:
                    new_dict[new_key] = {}  # Initialize nested dictionary
                new_dict[new_key][key] = value
            sample_hit_rank_2[sample] = new_dict

        sample_hit_rank_3 = {}
        # Keep top 3 values for each HLA group
        for sample, hlas in sample_hit_rank_2.items():
            sample_hit_rank_3[sample] = {}
            for hla, values in hlas.items():
                sample_hit_rank_3[sample][hla] = dict(sorted(values.items(), key=lambda x: x[1], reverse=True)[:top_n])

        # Initialize list to store results
        sample_hit_rank_low = []

        # Populate the list with sample, HLA type, rank, HLA, and probability
        for sample, hla_dict in sample_hit_rank_3.items():
            for hla_type, hla_values in hla_dict.items():
                sorted_hla = sorted(hla_values.items(), key=lambda x: x[1], reverse=True)
                for rank, (hla, probability) in enumerate(sorted_hla, 1):
                    sample_hit_rank_low.append([sample, hla_type.strip('*'), rank, hla, probability])

        # Create a DataFrame from the list
        top_hla_df = pd.DataFrame(
            sample_hit_rank_low, columns=['Sample', 'HLA_Type', 'Rank', 'HLA', 'Probability']
        )
        
        pred_dict = self.get_pred_dict(sample_hit_rank)
        # Add Final_decision column
        top_hla_df['Final_decision'] = 'False'
        for i in range(top_hla_df.shape[0]):
            sample = top_hla_df.loc[i, 'Sample']
            HLA = top_hla_df.loc[i, 'HLA']
            if HLA in pred_dict.get(sample, []):
                top_hla_df.loc[i, 'Final_decision'] = 'True'

        # Map additional metrics (AUC and F-score) to the DataFrame
        top_hla_df['AUC'] = top_hla_df['HLA'].map(self.hla_auc)
        top_hla_df['F_score'] = top_hla_df['HLA'].map(self.fscore_dict)

        return top_hla_df
