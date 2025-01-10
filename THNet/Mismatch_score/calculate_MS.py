#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os,sys
import pickle

class MS_calculation:
    def __init__(self):
        
        hla_list_path = os.path.join(os.path.dirname(__file__), 'parameter', 'hla_list.pkl')
        class1_distance_path = os.path.join(os.path.dirname(__file__), 'parameter', 'class1_distance.pkl')
        class2_distance_path = os.path.join(os.path.dirname(__file__), 'parameter', 'class2_distance.pkl')
        
        with open(hla_list_path, 'rb') as file:
            hla_list = pickle.load(file)
        with open(class1_distance_path, 'rb') as file:
            class1_distance = pickle.load(file)
        with open(class2_distance_path, 'rb') as file:
            class2_distance = pickle.load(file)
            
        self.hla_list = hla_list
        self.class1_distance = class1_distance
        self.class2_distance = class2_distance

    def _filter_hla(self, hla_data):
        return [hla for hla in hla_data if hla != 'X']

    def _is_valid_hla(self, hla_data):
        return all(hla in self.hla_list for hla in hla_data)

    def _calculate_mismatch_score(self, mismatch_pairs, distance_matrix, is_class1=True):
        mismatch_scores = {}
        for pair, mismatches in mismatch_pairs.items():
            mismatch_scores[pair] = 0
            tmp_dict = {}
            for mismatch in mismatches:
                hla1, hla2 = mismatch.split('_')
                if is_class1:
                    condition = len(hla1.split('*')[0]) == len(hla2.split('*')[0]) == 1
                else:
                    condition = hla1.split('*')[0] == hla2.split('*')[0] and len(hla2.split('*')[0]) == 4

                if condition:
                    n = distance_matrix[hla1][hla2]
                    if hla1 not in tmp_dict:
                        tmp_dict[hla1] = [n]
                    else:
                        tmp_dict[hla1].append(n)
            for hla, values in tmp_dict.items():
                mismatch_scores[pair] += sum(values) / len(values)
        return mismatch_scores

    def _get_mismatch_pairs(self, sample_pair):
        mismatch_pair = {}
        for pair, hlas in sample_pair.items():
            mismatch_pair[pair] = []
            for hla in sample_pair[pair]['donor']:
                if hla not in sample_pair[pair]['recipient']:
                    tmp_class1 = [hla + '_' + i for i in sample_pair[pair]['recipient'] if 'D' not in hla and 'D' not in i]
                    tmp_class2 = [hla + '_' + i for i in sample_pair[pair]['recipient'] if 'D' in hla and 'D' in i]
                    mismatch_pair[pair] += tmp_class1 + tmp_class2
            mismatch_pair[pair] = list(set(mismatch_pair[pair]))
        return mismatch_pair

    def process(self, input_df):
        sample_pair = {}
        sample_failed = []

        for row in range(input_df.shape[0]):
            sample = input_df.loc[row, 'TX_ID']
            rec_hla = self._filter_hla(input_df.loc[row, ['Rec_A_1', 'Rec_A_2', 'Rec_B_1', 'Rec_B_2',
                                                          'Rec_C_1', 'Rec_C_2', 'Rec_DQB1_1', 'Rec_DQB1_2',
                                                          'Rec_DRB1_1', 'Rec_DRB1_2']].tolist())
            don_hla = self._filter_hla(input_df.loc[row, ['Don_A_1', 'Don_A_2', 'Don_B_1', 'Don_B_2',
                                                          'Don_C_1', 'Don_C_2', 'Don_DQB1_1', 'Don_DQB1_2',
                                                          'Don_DRB1_1', 'Don_DRB1_2']].tolist())

            tmp_hla = rec_hla + don_hla
            if self._is_valid_hla(tmp_hla):
                sample_pair[sample] = {'recipient': rec_hla, 'donor': don_hla}
            else:
                sample_failed.append(sample)

        if sample_failed:
            sample_failed = [str(i) for i in sample_failed]
            print(f"Samples {', '.join(sample_failed)} contain HLA alleles that are not in the accepted HLA list.")

        input_df_filtered = input_df.loc[input_df['TX_ID'].isin(sample_pair.keys())]

        mismatch_pair = self._get_mismatch_pairs(sample_pair)

        mismatch_score_class1 = self._calculate_mismatch_score(mismatch_pair, self.class1_distance, is_class1=True)
        mismatch_score_class2 = self._calculate_mismatch_score(mismatch_pair, self.class2_distance, is_class1=False)

        output_df = pd.DataFrame({
            'TX_ID': list(mismatch_score_class1.keys()),
            'Class_I_MS': list(mismatch_score_class1.values()),
            'Class_II_MS': list(mismatch_score_class2.values())
        })

        return output_df
