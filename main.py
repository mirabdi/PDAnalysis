import argparse
import os
from pathlib import Path

from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np

from protein import Protein, AverageProtein
from deformation import Deformation


def parse_args():
    parser = argparse.ArgumentParser(
        prog="python main.py",
        description="Protein Deformation Analysis v0.0: software for analysing deformation between protein structures.")

    parser.add_argument("--protA", nargs="+", type=str, default='',
        help="Input Protein(s) A: Can be a path to a PDB, mmCIF, or coordinates file (.npy or .txt)." + \
             " Multiple paths can be provided separated by spaces (e.g., --protA path1 path2 path3).")

    parser.add_argument("--prot_listA", type=str, default='',
        help="Input Protein(s) A: Can be a path to a file multiple protein paths." + \
             "One protein path per line. Each protein path can be a PDB, mmCIF, or coordinates file (.npy or .txt).")

    parser.add_argument("--protB", nargs="+", type=str, default='',
        help="Input Protein(s) B: Can be a path to a PDB, mmCIF, or coordinates file (.npy or .txt)." + \
             " Multiple paths can be provided separated by spaces (e.g., --protA path1 path2 path3).")

    parser.add_argument("--prot_listB", type=str, default='',
        help="Input Protein(s) B: Can be a path to a file multiple protein paths." + \
             "One protein path per line. Each protein path can be a PDB, mmCIF, or coordinates file (.npy or .txt).")

    parser.add_argument("--max_bfactor", type=float, default=0.0,
        help="Maximum bfactor cutoff: exclude from calculations any atoms with B-factor higher than the cutoff." + \
             " Setting to zero turns this function off.")

    parser.add_argument("--min_plddt", type=float, default=0.0,
        help="Minimum pLDDT cutoff: exclude from calculations any atoms with pLDDT lower than the cutoff." + \
             " To avoid high deformation due to disordered residues, a value of 70 is suggested." + \
             " Setting to zero turns this function off.")
    
    parser.add_argument("--neigh_cut", type=float, default=13.0,
        help="Neighbor distance cutoff: distance in Angstroms within which residues are considered neighbors.")


    return parser.parse_args()


### Parses different input types, and return a list of paths
def parse_input_path_AB(prot, prot_list, lbl):
    if len(prot) and len(prot_list):
        raise Exception(f"ERROR! Please only specify either --prot{lbl} OR --prot_list{lbl}." + \
                         " Do not specify both at the same time")

    if len(prot):
        path_list = prot.copy()

    elif len(prot_list):
        if not os.path.exists(prot_list):
            raise FileNotFoundError(f"Path '{prot_list}' not found.")
        path_list = [l.strip('\n') for l in open(prot_list).readlines()]

    path_out = []
    for i, path in enumerate(path_list):
        if not os.path.exists(path):
            print(f"WARNING! Path {path} not found.\n\tContinuing without this path!")
        else:
            path_out.append(path)

    return path_out


def parse_input_paths(args):
    if len(args.protA) or len(args.prot_listA):
        pathA = parse_input_path_AB(args.protA, args.prot_listA, 'A')
    else:
        raise Exception("ERROR! No input files provided. Please specifiy at least --protA OR --prot_listA.")

    if len(args.protB) or len(args.prot_listB):
        pathB = parse_input_path_AB(args.protB, args.prot_listB, 'B')
    else:
        pathB = []

    if (not len(pathB)) & (len(pathA) == 1):
        raise Exception("ERROR! Only one protein file was provided. Please specify more than one file.")

    return pathA, pathB


def load_protein_object(path_list, **kwargs):
    if len(path_list) == 1:
        return Protein(path_list[0], **kwargs)

    elif len(path_list) > 1:
        return AverageProtein(path_list, **kwargs)


def load_protein_kwargs(args):
    kwargs_list = ["max_bfactor", "min_plddt", "neigh_cut"]
    return {k: getattr(args, k) for k in kwargs_list}


def main():
    args = parse_args()
    pathA, pathB = parse_input_paths(args)
    protein_kwargs = load_protein_kwargs(args)

    if len(pathB):
        protA = load_protein_object(pathA, **protein_kwargs)
        protB = load_protein_object(pathB, **protein_kwargs)
        deform = Deformation(protA, protB)
        deform.calculate_deformation()
        np.savetxt("strain.txt", deform.strain)

    else:
        protA = load_protein_object(pathA, **protein_kwargs)
        


if __name__ == "__main__":
    main()
    


