#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from THNet.HLA_inference.model_prediction import Model_prediction
from THNet.HLA_inference.HLA_inference import HLA_inference
from THNet.Mismatch_score.calculate_MS import MS_calculation
import pandas as pd, numpy as np
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

"""Code to execute the functionality of THNet

    Copyright (C) 2024 Mingyao Pan and Bo Li

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""

def main():
    parser = argparse.ArgumentParser(description="Initialize THNet")

    # Sub-command parser
    subparsers = parser.add_subparsers(dest="command", required=True, help="Select function (HLA_inference/Mismatch_score)")

    # HLA_inference sub-command
    parser_a = subparsers.add_parser("HLA_inference", help="Execute HLA_inference")
    parser_a.add_argument('-i', '--input_file', dest='input_path', metavar='PATH/TO/FILE',
                          help='Load input file that contains CDR3 beta sequences and V gene families from PATH/TO/FILE', required=True)
    parser_a.add_argument('-o', '--output_file', dest='output_folder', metavar='PATH/TO/FOLDER',
                          help='Write model output to PATH/TO/FOLDER', required=True)
    parser_a.add_argument('-n', '--Top_HLA_n', default=3, dest='top_n', metavar='n',
                          choices=[1, 2, 3, 4, 5], help='Output the top n most probable HLA alleles for each HLA type')

    # Mismatch_score sub-command
    parser_b = subparsers.add_parser("Mismatch_score", help="Execute Mismatch_score")
    parser_b.add_argument('-i', '--input_file', dest='input_path', metavar='PATH/TO/FILE',
                          help='Load input file that contains HLA alleles of both donor and recipient from PATH/TO/FILE', required=True)
    parser_b.add_argument('-o', '--output_file', dest='output_folder', metavar='PATH/TO/FOLDER',
                          help='Write model output to PATH/TO/FOLDER', required=True)

    # Parse command-line arguments
    args = parser.parse_args()

    # Check if input file exists
    if not os.path.exists(args.input_path):
        raise FileNotFoundError(f"Input file not found: {args.input_path}")

    # Ensure the output directory exists
    os.makedirs(args.output_folder, exist_ok=True)

    # Read input file
    try:
        input_file = pd.read_csv(args.input_path)
    except Exception as e:
        raise ValueError(f"Failed to load input file: {args.input_path}. Error: {e}")

    if args.command == "HLA_inference":
        # Perform HLA inference
        model_prediction = Model_prediction()
        sample_hit_rank = model_prediction.Get_prediction(input_file)

        processor = HLA_inference()

        # Write HLA_inference.csv
        pred_result = processor.hla_inference_df(sample_hit_rank)
        pred_result_path = os.path.join(args.output_folder, 'HLA_inference.csv')
        pred_result.to_csv(pred_result_path, index=False)

        # Write Top_hlas.csv
        top_hla_df = processor.create_top_hla_df(sample_hit_rank, top_n=args.top_n)
        top_hla_path = os.path.join(args.output_folder, 'Top_hlas.csv')
        top_hla_df.to_csv(top_hla_path, index=False)

    elif args.command == "Mismatch_score":
        # Perform Mismatch Score calculation
        ms_calculation = MS_calculation()
        mismatch_pair = ms_calculation.process(input_file)

        # Write TX_Mismatch_score.csv
        mismatch_path = os.path.join(args.output_folder, 'TX_Mismatch_score.csv')
        mismatch_pair.to_csv(mismatch_path, index=False)


if __name__ == '__main__':
    main()
