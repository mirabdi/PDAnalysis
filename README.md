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
* pandas
* Biopython

### Installing

* pip install PDAnalysis

### Executing program

```
pdanalysis --protA [path(s) to protein A] --protB [path(s) to protein B]
```

## Help

Any advise for common problems or issues.
```
command to run if program contains helper info
```

## Authors

[John M. McBride](https://github.com/jomimc)
[Amirbek Abdirasulov](https://github.com/amirbek)

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

