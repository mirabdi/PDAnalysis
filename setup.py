#!/usr/bin/env python
# PDAnalysis (Protein Deformation Analysis)
# Copyright (C) 2023 John Michael McBride
#
# Released under MIT License
#
# Please cite use of PDAnalysis in published work:
#
# J. M. McBride, K. Polev, A. Abdirasulov, V. Reinharz, B. A. Grzybowski,
# T. Tlusty.
# AlphaFold2 can predict single mutation effects, PRL (2023)
#

"""Setuptools-based setup script for PDAnalysis.

See ./requirements.txt for a list of dependencies.

For a basic installation just type the command::

  python setup.py install

"""

from setuptools import setup, find_packages
from pathlib import Path
import sys


# Check Python version
if sys.version_info[:2] < (3, 7):
    print('PDAnalysis requires Python 3.7+. Python {0:d}.{1:d} detected'.format(*
          sys.version_info[:2]))
    print('Please upgrade your version of Python.')
    sys.exit(-1)


if __name__ == "__main__":

    HERE = Path(__file__).parent.resolve()
    LONG_DESCRIPTION = (HERE / "README.md").read_text(encoding="utf-8")

    CLASSIFIERS = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: C',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]

    INSTALL_REQUIRES = [        
        'numpy>=1.22.3',
        'biopython>=1.80',
        'pandas>=1.4.4',
        'scipy>=1.5.0'
    ]

    ENTRY_POINTS = {
        "console_scripts": [
            "pdanalysis=PDAnalysis.main:main"
        ]
    }

    setup(name='PDAnalysis',
          version="0.0.1",
          description='Software for analysing deformation between protein structures.',
          long_description=LONG_DESCRIPTION,
          long_description_content_type='text/markdown',
          author='John M. McBride',
          author_email='jmmcbride@protonmail.com',
          url='https://www.mdanalysis.org',
          classifiers=CLASSIFIERS,
          provides=['PDAnalysis'],
          packages=find_packages(),
          python_requires='>=3.7',
          install_requires=INSTALL_REQUIRES,
          keywords="science biology biophysics molecular structural",
          entry_points=ENTRY_POINTS
        )



