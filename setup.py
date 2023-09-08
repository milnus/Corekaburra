#!/usr/bin/env python

from distutils.core import setup

LONG_DESCRIPTION = \
'''Corekaburra looks at the gene synteny across genomes used to build a pan-genome. Using syntenic information Corekaburra 
identifies regions between core gene clusters. Regions are described in terms of their content of accessory gene clusters 
and distance between core genes. Information from neighboring core genes is further used to identify stretches of core  
gene clusters throughout the pan-genome that appear in all genomes given as input. Corekaburra is compatible with outputs 
from standard pan-genome pipelines: Roary and Panaroo.'''


setup(
    name='Corekaburra',
    version='0.0.5',
    author='Magnus Ganer Jespersen',
    author_email='magnus.ganer.j@gmail.com',
    packages=['Corekaburra'],
    package_dir={'Corekaburra': 'Corekaburra'},
    entry_points={
        'console_scripts': ['Corekaburra = Corekaburra.__main__:main']
    },
    url='https://github.com/milnus/Corekaburra',
    license='LICENSE',
    description=('A commandline bioinformatics tool to utilize syntenic information from genomes in the context of pan-genomes'),
    long_description=(LONG_DESCRIPTION),
    install_requires=["biopython==1.79", "networkx>=2.6.3", "gffutils>=0.10.1", "numpy>=1.23.4"],
    keywords=['Genomics', 'pan-genome', 'bacteria', 'prokaryotes', 'bioinformatics'],
    classifiers=[
        'Programming Language :: Python :: 3.9',
        'License :: OSI Approved :: MIT License',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Development Status :: 4 - Beta']
)
