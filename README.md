# Protein Deformation Analysis (PDAnalysis)

Python package for calculating deformation between protein structures.

## Description

Calculates local similarity:
* Local Distance Difference Test (LDDT)

Calculates local deformation:
* Local Distance Difference (LDD)
* Neighborhood Distance
* Effective Strain
* Shear Strain
* Non-affine Strain

Calculates global deformation:
* Root-mean-squared-deviation (RMSD)


## Getting Started

### Dependencies

* python 3.7<=
* numpy
* scipy
* pandas
* Biopython


### Installing

* Download this repository
* pip install PDAnalysis (in development - not yet ready)

### Executing program

```
python main.py --protA [path(s) to protein A] --protB [path(s) to protein B]
```

Protein Input Types:
* "--protA"         :: Path to a single protein, or multiple paths separated by spaces. Allowed input types: .pdb, .cif, .npy, .txt.
* "--protB"         :: Path to a single protein, or multiple paths separated by spaces. Allowed input types: .pdb, .cif, .npy, .txt.
* "--prot\_listA"   :: Path to a simple text file with a list of paths to proteins.
* "--prot\_listB"   :: Path to a simple text file with a list of paths to proteins.

Optional Arguments:
* "--method"        :: Deformation method to use. Default is "strain". "all" calculates everything. Multiple methods can be used if separated by spaces. Allowed methods: "mut\_dist", "strain", "shear", "non-affine", "ldd", "lddt", "rmsd". Note that "rmsd" only works for comparisons between single structures (no averaging); running "rmsd" with multple input structures will only use one protein from A and one protein from B.
* "--min\_plddt"    :: Minimum pLDDT cutoff: exclude from calculations any atoms with pLDDT lower than the cutoff. Recommended value = 70.
* "--max\_bfactor"  :: Maximum bfactor cutoff: exclude from calculations any atoms with B-factor higher than the cutoff.
* "--neigh\_cut"    :: Neighbor cutoff distance in Angstroms.
* "--lddt\_cutoffs" :: Specify the distance cutoffs used in LDDT calculation. Default = (0.5, 1, 2, 4).


Examples:  
Calculates "Effective Strain" between two PDB conformations, and ignore residues with pLDDT < 70:
```
python main.py --protA test_data/Lysozyme/AF-P61626-F1-model_v4.pdb --protB test_data/Lysozyme/AF-P79180-F1-model_v4.pdb --min_plddt 70
```

Calculates "LDDT" using one PDB file and two PDB files (an averaged configuration) and use a neighbor cutoff distance of 10 Angstroms.
```
python main.py --protA test_data/GFP/GFP_mutant.pdb --protB test_data/GFP/GFP_WT.pdb test_data/GFP/GFP_WT_ver2.pdb --neigh_cut 10 --method lddt
```

Calculates "Effective Strain" and "Shear Strain" using two lists of PDB files (averaging over all proteins in each list), use a neighbor cutoff distance of 12 Angstroms, and ignore residues with pLDDT < 70.
```
python main.py --protA test_data/GFP/mut_all.txt --protB test_data/GFP/wt_all.txt --neigh_cut 12 --min_plddt 70 --method strain shear
```

Calculates all measures between two PDB conformations, and uses alternative LDDT cutoff values:
```
python main.py --protA test_data/Lysozyme/AF-P61626-F1-model_v4.pdb --protB test_data/Lysozyme/AF-P79180-F1-model_v4.pdb --lddt_cutoffs 0.125 0.25 0.5 1
```

## Help

For more detailed help on running code:
```
python main.py --help
```

## Authors

* [John M. McBride](https://github.com/jomimc)
* [Amirbek Abdirasulov](https://github.com/amirbek)
* [Tsvi Tlusty](http://www.sns.ias.edu/~tlusty/index.html)

## Version History

* 0.0
    * Pre-Release

## License

This project is licensed under the MIT License - see the LICENSE.md file for details

## Acknowledgments

Inspiration, code snippets, etc.
* LDDT
* Shear Strain (Jacques Rougemont, paper)
* Non-affine Strain (paper, matscipy)

