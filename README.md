# THNet (TCR-based HLA similarity mapping network)
THNet is a Python 3.11 software designed to infer HLA haplotypes from T-cell beta chain repertoire datasets that can calculate mismatch scores (MS) by taking the HLA allele compositions of both donor and recipient as input to predict the transplantation outcome. The model can infer 208 HLA alleles based on the T-cell beta chain repertoire, with higher accuracy for common alleles. 

THNet is developed and maintained by [Li lab at the University of Pennsylvania](https://lilab-utsw.org). Please direct your questions regarding THNet to Mingyao Pan: mingyaop@seas.upenn.edu.

THNet is written in Python3.11, with the following dependencies:

- [Panda]
- [numpy]
- [tqdm]
- [scikit-learn]

## Installation
THNet is available on PyPI and can be downloaded and installed through pip:

 ```pip install THNet```

THNet is also available on [GitHub](https://github.com/Mia-yao/THNet/tree/main). The command line entry points can be installed by using the setup.py script:

 ```$python setup.py install ```

 Directory architecture:
```
 THNet/
├── README.md
├── LICENSE
├── setup.py
├── MANIFEST.in
├── THNet/                 
│   ├── __init__.py
│   ├── load_model.py
│   ├── HLA_inference/   
│   │   ├── __init__.py
│   │   ├── HLA_inference.py
│   │   ├── model_prediction.py
│   │   ├── example/
│   │   │   └── input_example.csv
│   │   ├── models/
│   │   │   ├── models_1.pkl  
│   │   │   └── models_2.pkl 
│   │   └── parameter/
│   │       ├── fscore_dict.pkl
│   │       ├── hla_auc.pkl
│   │       ├── hla_list.pkl
│   │       ├── hla_threshold.pkl      
│   │       └── v_gene_list.pkl
│   └── Mismatch_score/    
│       ├── __init__.py
│       ├── calculate_MS.py
│       ├── example/
│       │   └── input_example.csv
│       └── parameter/
│           ├── class1_distance.pkl  
│           ├── class2_distance.pkl 
│           └── hla_list.pkl
```

## Usage

Type `THNet --help` to display all the command line options:

THNet has two functions: `HLA_inference` and `Mismatch_score`. 

### HLA_inference
`HLA_inference` is for the inference of HLA based on one's T cell beta chain repertoire

|Commands|Description|
|--|--|
|`-i, --input_file`|Load input file that contains CDR3 beta sequences and V gene families from PATH/TO/FILE.| 
|`-o, --output_file`|Write model output to PATH/TO/FOLDER|
|`-n --Top_HLA_n`|Output the top n most probable HLA alleles for each HLA type. Default 3. The valid value of n is 1 to 5|

* Input data format
  
The input file of HLA_inference is a .csv file (separated by delimiter ',') containing three columns: sample, cdr3, and v_gene. Note that the format of the V gene has to be: TRBVXX-XX (IMGT format). The vaild v gene list can be checked at `THNet/HLA_inference/parameter/v_gene_list.pkl`
```
cdr3,v_gene,sample
CAWSRGGVTGELFF,TRBV30,Sample1
CASKPMVNEQFF,TRBV19,Sample1
CASSLGAGLQETQYF,TRBV13,Sample1
CASSLSSGSSYNEQFF,TRBV27,Sample1
CASNAGLRDTQYF,TRBV2,Sample1
CASSAGTVVGNTIYF,TRBV5-1,Sample1
```
An example of the input files for Mismatch_score can be referred to at `THNet/HLA_inference/example/input_example.csv`

Note: A single sample should have around 10,000 TCR sequences for better model performance. Samples with very few TCRs will yield minimal HLA inference results. 

* Output data format
  
The output includes two files: `HLA_inference.csv`, which contains the final HLA predictions, and `Top_hlas.csv`, listing the top n HLA alleles with the highest probabilities for each HLA type.

* Demo usage
  
`THNet HLA_inference -i input_file_path/input.csv -o output_folder_path -n 4`

This command line takes `input.csv` as input data and outputs the result files `HLA_inference.csv` and `Top_hlas.csv` in the output_folder_path folder.

Note: The input data is provided in a file, and the output_folder_path specifies where the output files will be stored, as there are multiple output files.

### Mismatch_score
`Mismatch_score` is for the calculation of the mismatch scores (MS) by taking the HLA allele compositions of both donor and recipient as input to predict the outcome of transplantation.

|Commands|Description|
|--|--|
|`-i, --input_file`|Load input file that contains CDR3 beta sequences and V gene families from PATH/TO/FILE.| 
|`-o, --output_file`|Write model output to PATH/TO/FOLDER|

* Input data format
The input file of HLA_inference is a .csv (separated by delimiter ',') file containing 17 columns: TX_ID,Rec_A_1,Rec_A_2,Rec_B_1,Rec_B_2,Rec_C_1,Rec_C_2,Rec_DQB1_1,Rec_DQB1_2,Rec_DRB1_1,Rec_DRB1_2,
Don_A_1,Don_A_2,Don_B_1,Don_B_2,Don_C_1,Don_C_2,Don_DQB1_1,Don_DQB1_2,Don_DRB1_1,Don_DRB1_2。
An example of the input files for Mismatch_score can be referred to at `THNet/Mismatch_score/example/input_example.csv`

Note: Some HLA information may be missing, but a capitalized 'X' should be placed in the corresponding position of the table as a placeholder. All 17 columns are required. The mismatch score will be calculated using the available HLA alleles for that transplantation pair. If all class I or class II HLA alleles of the donor or recipient are missing, the `Class_I_MS` or `Class_II_MS` will be set to zero.

Note: So far, there are 208 valid HLA alleles (including both class I and II). The valid HLA alleles can be checked in the file `THNet/Mismatch_score/parameter/hla_list.pkl`. Any HLA allele not present in this list will trigger an error message.

* Output data format
  
The output includes one file: `TX_Mismatch_score.csv`, which contains three columns: `TX_ID`,	`Class_I_MS`, and `Class_II_MS`. `Class_I_MS` and `Class_II_MS` indicate the mismatch score for HLA class I and HLA class II of the given TX_ID respectively.

* Demo usage
  
`THNet Mismatch_score -i input_file_path/input.csv -o output_folder_path`

This command line takes `input.csv` as input data and outputs the result files `TX_Mismatch_score.csv` in the output_folder_path folder.
