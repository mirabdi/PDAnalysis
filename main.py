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

    ### Protein input files

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


    ### Parameters controlling which residues to include in deformation calculations

    parser.add_argument("--max_bfactor", type=float, default=0.0,
        help="Maximum bfactor cutoff: exclude from calculations any atoms with B-factor higher than the cutoff." + \
             " Setting to zero turns this function off.")

    parser.add_argument("--min_plddt", type=float, default=0.0,
        help="Minimum pLDDT cutoff: exclude from calculations any atoms with pLDDT lower than the cutoff." + \
             " To avoid high deformation due to disordered residues, a value of 70 is suggested." + \
             " Setting to zero turns this function off.")
    
    parser.add_argument("--neigh_cut", type=float, default=13.0,
        help="Neighbor distance cutoff: distance in Angstroms within which residues are considered neighbors.")
    

    ### Other

    parser.add_argument("-m", "--method", nargs="+", type=str, default='',
        help="Deformation method: 'all' will automatically calculate everything. Multiple methods may be provided separated by spaces. Accepted methods:" + \
             "\n\t'mut_dist' :: Distance from nearest mutated residue" + \
             "\n\t'strain' :: Effective Strain" + \
             "\n\t'shear' :: Shear Strain" + \
             "\n\t'non-affine' :: Non-Affine Strain" + \
             "\n\t'ldd' :: Local Distance Difference (LDD)" + \
             "\n\t'lddt' :: Local Distance Difference Test (LDDT)" + \
             "\n\t'neighborhood_dist' :: Neighborhood Distance" + \
             "\n\t'rmsd' :: Root-mean-squared-deviation (RMSD)")

    parser.add_argument("--lddt_cutoffs", type=float, nargs="+", default=[0.5, 1, 2, 4],
        help="Neighbor distance cutoff: distance in Angstroms within which residues are considered neighbors.")

    parser.add_argument("-v", "--verbose", default=False, action='store_true',
        help="Print out information in addition to WARNINGS and ERRORS.")

    parser.add_argument("-o", "--output", default='output.csv', type=str,
        help="Path to output file.")


    return parser.parse_args()



### Parses different input types, and return a list of paths
def parse_input_path_AB(prot, prot_list, lbl):
    if len(prot) and len(prot_list):
        raise Exception(f"ERROR! Please only specify either --prot{lbl} OR --prot_list{lbl}." + \
                         " Do not specify both at the same time")

    # Get list of paths
    if len(prot):
        path_list = prot.copy()

    elif len(prot_list):
        if not os.path.exists(prot_list):
            raise FileNotFoundError(f"Path '{prot_list}' not found.")
        path_list = [l.strip('\n') for l in open(prot_list).readlines()]

    # Check that paths exist
    path_out = []
    for i, path in enumerate(path_list):
        if not os.path.exists(path):
            print(f"WARNING! Path {path} not found.\n\tContinuing without this path!")
        else:
            path_out.append(path)

    return path_out


### Parse input to get lists of paths
def parse_input_paths(args):
    # Do not allow path lists to be defined using both input types simulaneously
    if len(args.protA) or len(args.prot_listA):
        pathA = parse_input_path_AB(args.protA, args.prot_listA, 'A')
    else:
        raise Exception("ERROR! No input files provided. Please specifiy at least --protA OR --prot_listA.")

    # It is permitted to not pass any arguments for protB
    if len(args.protB) or len(args.prot_listB):
        pathB = parse_input_path_AB(args.protB, args.prot_listB, 'B')
    else:
        pathB = []

    # More than one protein file must be submitted
    if (not len(pathB)) & (len(pathA) == 1):
        raise Exception("ERROR! Only one protein file was provided. Please specify more than one file.")

    return pathA, pathB


### Load a Protein if the input is only one configuration,
### load an AverageProtein if the input consists of multiple configurations
def load_protein_object(path_list, **kwargs):
    if len(path_list) == 1:
        return Protein(path_list[0], **kwargs)

    elif len(path_list) > 1:
        return AverageProtein(path_list, **kwargs)


def load_protein_kwargs(args):
    kwargs_list = ["max_bfactor", "min_plddt", "neigh_cut"]
    return {k: getattr(args, k) for k in kwargs_list}


def load_deformation_kwargs(args):
    kwargs_list = ["method", "max_bfactor", "min_plddt", "neigh_cut", "verbose", "lddt_cutoffs"]
    return {k: getattr(args, k) for k in kwargs_list}


def main():
    args = parse_args()
    pathA, pathB = parse_input_paths(args)
    protein_kwargs = load_protein_kwargs(args)
    deform_kwargs = load_deformation_kwargs(args)

    if len(pathB):
        protA = load_protein_object(pathA, **protein_kwargs)
        protB = load_protein_object(pathB, **protein_kwargs)
        deform = Deformation(protA, protB, **deform_kwargs)
        deform.calculate_deformation()
        deform.save_output(args.output)

    else:
        protA = load_protein_object(pathA, **protein_kwargs)
        


if __name__ == "__main__":
    main()
    


