from setuptools import setup, find_packages
import sys

sys.path.append("./THNet")

def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name='THNet',
    version='1.0.2',
    description='HLA typing based on T-cell beta chain repertoires and HLA mismatch score calculation.',
    long_description=readme(),
    long_description_content_type='text/markdown', 
    url='https://github.com/Mia-yao/THNet',
    author='Mingyao Pan',
    author_email='mingyaop@seas.upenn.edu',
    license='MIT License',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Healthcare Industry',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.11',
    ],
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'tqdm',
        'scikit-learn'
    ],
    package_data={
        'THNet': [
            'HLA_inference/example/*',
            'HLA_inference/models/*',
            'HLA_inference/parameter/*',
            'Mismatch_score/example/*',
            'Mismatch_score/parameter/*',
        ],
    },
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'THNet=THNet.Load_model:main',  
        ],
    },
    zip_safe=False
)
