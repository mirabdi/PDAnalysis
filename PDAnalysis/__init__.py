"""Protein Deformation Analysis (PDAnalysis)

This package reads coordinates of two proteins (or two sets of proteins) A and B,
and calculates deformation between them.

"""

from .protein import Protein, AverageProtein
from .deformation import Deformation
from . import pdb_parser, utils



