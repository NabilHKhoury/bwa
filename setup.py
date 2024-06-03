from setuptools import setup, find_packages

setup(
    name='bwalign',
    version='1.6.6',    
    description='A burrows-wheeler seed and extend aligner',
    url='https://github.com/NabilHKhoury/bwalign',
    author=['Adrian Layer, Nabil Khoury, Yasmin Jabir'],
    license='MIT',
    packages=find_packages(),
    install_requires=['biopython',
                      'pysam',
                      'tqdm',                
                      ],
    entry_points={
        'console_scripts': [
            'bwalign=bwalign.main:main',
        ],
    },

    classifiers=[
        'Intended Audience :: Science/Research',  
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
)