from setuptools import setup, find_packages
import sys
sys.path.append("./THNet")

def readme():
    with open('README.md') as f:
        return f.read()

data_files_to_include = [('', ['README.md', 'LICENSE'])]

setup(
    name='THNet',
    version='1.0.0',
    description='HLA typing based on T-cell beta chain repertoires and HLA mismatch score calculation.',
    long_description='THNet is a Python 3.11 software designed to infer HLA haplotypes from T-cell beta chain repertoire datasets. '
                     'It accepts a .csv input file containing three columns: sample, cdr3, and v_gene. The model can infer 208 HLA alleles based on the T-cell beta chain repertoire, '
                     'with higher accuracy for common alleles. The output includes two files: HLA_inference.csv, which contains the final HLA predictions, and Top_hlas.csv, '
                     'listing the top n HLA alleles with the highest probabilities for each HLA type. Additionally, THNet can calculate mismatch scores (MS) by taking the HLA allele '
                     'compositions of both donor and recipient as input, outputting HLA class I and class II mismatch scores based on the HLA distances computed by the model.',
    long_description_content_type='text/plain', 
    url='https://github.com/Mia-yao/THNet',
    author='Mingyao Pan',
    author_email='mingyaop@seas.upenn.edu',
    license='GPLv3',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.11',
    ],
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'tqdm',
    ],
    package_data={
        'THNet': ['HLA_inference/example/*', 'Mismatch_score/example/*'],  
    },
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'THNet=THNet.Load_model:main',
        ],
    },
    zip_safe=False
)
