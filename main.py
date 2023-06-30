import argparse
from pathlib import Path

from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np

from protein import Protein, AverageProtein
from deformation import Deformation


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=str)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    print(args.inputs)
    


