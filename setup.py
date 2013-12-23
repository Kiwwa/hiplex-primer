#!/usr/bin/env python

from distutils.core import setup

setup(
    name='hiplex-primer',
    version='0.1.0',
    author='Bernie Pope',
    author_email='bjpope@unimelb.edu.au',
    packages=['Primer_design'],
    scripts=['Primer_design/Primer_design.py'],
    entry_points={
        'console_scripts': ['hiplex-primer = Primer_design:main']
    },
    url='https://github.com/bjpop/hiplex-primer',
    licence='LICENSE.txt',
    description=(
        'Hiplex-Primer: A tool for generating primers for the Hi-Plex \
         targeted, multiplexed DNA sequencing strategy.'),
    install_requires=['biopython >= 1.62'],
)