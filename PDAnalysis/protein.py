from collections import Counter
from pathlib import Path

from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist

from pdb_parser import parse_pdb_coordinates, parse_mmcif_coordinates, load_and_fix_pdb_data
from utils import rotate_points


class Protein:
    """ Protein object description:

    The Protein object accepts a protein structure as input (PDB file, CIF file, atomic coordinates),
    and represents the protein conformation as a list of local neighborhoods per residue.

    Attributes
    ----------
    seq_len : int
        length of the protein sequence

    sequence : np.ndarray(seq_len)
        single-letter amino-acid sequence

    coord : np.ndarray(seq_len, 3)
        xyz coordinates of alpha-carbons, with nan values for missing residues
        and for residues excluded based on bfactor/pLDDT
        
    coord_raw : np.ndarray(seq_len, 3)
        xyz coordinates of alpha-carbons, exactly as read from file

    idx : np.ndarray(seq_len)
        residue indices

    bfactor : np.ndarray(seq_len)
        Bfactor (temperature factor); mirrors plddt

    plddt : np.ndarray(seq_len)
        pLDDT (AlphaFold-predicted confidence score); mirrors bfactor

    dist_mat : np.ndarray(seq_len, seq_len)
        alpha-carbon distance matrix

    neigh_idx : list
        list of neighbor indices for each residue

    neigh_tensor: list
        list of neighborhood tensors for each residue


    Methods
    -------

    get_local_neighborhood
        Gets the indices of neighbors for each residue (neigh_idx),
        and calculates neighborhood tensors (neigh_tensor).

    """
    def __init__(self, inp_obj, **kwargs):
        """
        args
        ----------
        inp_obj : (str, Path, np.ndarray)
            > path (str, Path) to file:
                  allowed filetypes:
                      PDB / CIF file (.pdb, .ent, .cif)
                      coordinate text file (.txt, .dat)
                      coordinate numpy binary file (.npy)
            > coordinates (np.ndarray)

        kwargs
        ----------

        max_bfactor : float : default = 0.0
            if nonzero, excludes residues with zscore(bfactor) greater than max_bfactor

        min_plddt : float : default = 0.0
            if nonzero, excludes residues with pLDDT lower than min_plddt

        neigh_cut : float : default = 13.0
            cutoff radius used to determine neighborhood

        chain : str : default = ''
            if not '', only parses coordinates for a specific chain

        fix_pdb : bool : default = False
            if True, uses the SEQRES data to assign residue indices;
            neighborhoods of missing residues are stored as empty lists;
            this facilitates comparison of PDB structures that are missing different sets of residues

        delimiter : str : default = ' '
            delimiter for loading coordinates in text files;
            passed to np.loadtxt

        skip_rows : int : default = 0
            number of rows to skip at the top of text files;
            passed to np.loadtxt
            
        """
        # Parameters
        self.max_bfactor = kwargs.get('max_bfactor', 0.0)
        self.min_plddt = kwargs.get('min_plddt', 0.0)
        self.neigh_cut = kwargs.get('neigh_cut', 13.0)
        self.chain = kwargs.get('chain', '')
        self.fix_pdb =  kwargs.get('fix_pdb', False)
        self.delimiter = kwargs.get('delimiter', ' ')
        self.skip_rows = kwargs.get('skip_rows', 0)

        # Defaults
        self.sequence = None
        self.idx = None
        self.plddt = None
        self.bfactor = None
        self.dist_mat = None
        self.neigh_idx = []
        self.neigh_tensor = []
        self.seq_len = None         # Number of residues

        # Load data
        self._parse_input(inp_obj)
        self.get_local_neighborhood()


    def _parse_input(self, inp_obj):
        """
        Parses a range of input object types (.pdb, .ent, .cif, .txt, .dat, .npy)

        Extracts: coordinates, residue indices, amino-acid sequence, bfactor/pLDDT
        """

        if isinstance(inp_obj, (str, Path)):
            ext = Path(inp_obj).suffix

            ### Read PDB / mmCIF file
            if (ext in ['.pdb', '.cif', '.ent']):
                self._load_data_from_path(inp_obj, ext)

            ### Read text coordinate file
            elif (ext in ['.txt', '.dat']):
                try:
                    print(f"WARNING! Ambiguity in file format {ext}. Attempting to read coordinates.")
                    self.coord = np.loadtxt(inp_obj, delimiter=self.delimiter, skip_rows=self.skip_rows)
                    self.idx = np.arange(len(self.coord))
                except:
                    raise Exception("Could not read coordinates from input file!")

            ### Read nump binary coordinate file
            elif (ext == '.npy'):
                self.coord = np.load(inp_obj)
                self.idx = np.arange(len(self.coord))

        elif isinstance(inp_obj, np.ndarray):
            self.coord = inp_obj
            self.idx = np.arange(len(self.coord))

        else:
            raise Exception(f"Input object type {type(inp_obj)} is not supported.")

        if not isinstance(self.sequence, type(None)):
            self.seq_len = len(self.sequence)
        else:
            self.seq_len = len(self.coord)


    def _load_data_from_path(self, path, ext):
        """Parses data from a path, depending on the file type"""
        if self.fix_pdb:
            data = load_and_fix_pdb_data(path, self.chain)
            # (indices in PDB files are inconsistent, so we reorder
            #  the indices so that they start from zero)

        else:
            if ext in [".pdb", ".ent"]:
                data = parse_pdb_coordinates(path, self.chain)
            elif ext in [".cif"]:
                data = parse_mmcif_coordinates(path, self.chain)

            # Reorder indices so that they start from zero
            # (by default, AF-predicted structures start counting at one)
            if data[1][0] != 0:
                data[1] = data[1] - data[1][0]

        self.coord_raw = data[0]
        self.idx = data[1]
        self.sequence = data[2]
        self.bfactor = data[3]
        self.plddt = data[3].copy()

        # Fill in 'disordered' (high bfactor, low plddt) with NaN values
        self._update_nan_coords()


    # Fill in 'disordered' (high bfactor, low plddt) with NaN values
    def _update_nan_coords(self):
        """Constructs a coordinate array that includes nan values for
        missing residues and for residues excluded based on bfactor/pLDDT"""

        self.coord = np.zeros((len(self.sequence), 3)) * np.nan
        self.coord[self.idx] = self.coord_raw.copy()

        bfactor = np.zeros(len(self.sequence)) * np.nan
        bfactor[self.idx] = self.bfactor.copy()
        self.bfactor = bfactor

        if self.max_bfactor > 0.0:
            bfactor_zscore = (self.bfactor - np.nanmean(self.bfactor)) / np.nanstd(self.bfactor)
            self.coord[(bfactor_zscore > self.max_bfactor) | (np.isnan(self.bfactor))] = np.nan
        elif self.min_plddt > 0.0:
            self.coord[self.plddt <= self.min_plddt] = np.nan


    def _get_dist_mat(self):
        """Calculates alpha-carbon distance matrix"""
        self.dist_mat = cdist(self.coord, self.coord)
    

    def get_local_neighborhood(self):
        """Extracts local neighborhoods for each residue"""
        if not isinstance(self.dist_mat, np.ndarray):
            self._get_dist_mat()
        self.neigh_idx = [np.where((d > 0) & (d <= self.neigh_cut) & (np.isfinite(d)))[0] for d in self.dist_mat]
        self._calculate_neighbor_tensor()


    def _calculate_neighbor_tensor(self):
        """Calculates the neighborhood tensor for each residue"""
        self.neigh_tensor = []
        for i, idx in enumerate(self.neigh_idx):
            self.neigh_tensor.append(self.coord[idx] - self.coord[[i]])



class AverageProtein:
    """ AverageProtein object description:

    The AverageProtein object accepts a protein structure as input (PDB file, CIF file, atomic coordinates),
    and represents the protein conformation as a list of local neighborhoods per residue.

    Attributes
    ----------
    seq_len : int
        length of the protein sequence

    num_repeat : int
        number of proteins used to create an AverageProtein

    sequence : np.ndarray(seq_len)
        single-letter amino-acid sequence

    coord : np.ndarray(seq_len, 3)
        xyz coordinates of alpha-carbons, with nan values for missing residues
        and for residues excluded based on bfactor/pLDDT
        
    coord_raw : np.ndarray(seq_len, 3)
        xyz coordinates of alpha-carbons, exactly as read from file

    bfactor : np.ndarray(seq_len)
        Bfactor (temperature factor); mirrors plddt

    plddt : np.ndarray(seq_len)
        pLDDT (AlphaFold-predicted confidence score); mirrors bfactor

    dist_mat : np.ndarray(seq_len, seq_len)
        alpha-carbon distance matrix

    neigh_idx : list
        list of neighbor indices for each residue

    neigh_tensor: list
        list of average neighborhood tensors for each residue


    Methods
    -------

    get_average_structure
        Gets the indices of neighbors for each residue (neigh_idx),
        and calculates neighborhood tensors (neigh_tensor).

    recalculate_average_structure
        Recalculates the average structure after changing parameters.

    """
    def __init__(self, inp_obj, **kwargs):
        """
        args
        ----------
        inp_obj : list of (Protein, str, Path, np.ndarray)
            > Protein objects
            > path (str, Path) to file:
                  allowed filetypes:
                      PDB / CIF file (.pdb, .ent, .cif)
                      coordinate text file (.txt, .dat)
                      coordinate numpy binary file (.npy)
            > coordinates (np.ndarray)

        kwargs
        ----------

        max_bfactor : float : default = 0.0
            if nonzero, excludes residues with zscore(bfactor) greater than max_bfactor

        min_plddt : float : default = 0.0
            if nonzero, excludes residues with pLDDT lower than min_plddt

        neigh_cut : float : default = 13.0
            cutoff radius used to determine neighborhood

        chain : str : default = ''
            if not '', only parses coordinates for a specific chain

        fix_pdb : bool : default = False
            if True, uses the SEQRES data to assign residue indices;
            neighborhoods of missing residues are stored as empty lists;
            this facilitates comparison of PDB structures that are missing different sets of residues

        delimiter : str : default = ' '
            delimiter for loading coordinates in text files;
            passed to np.loadtxt

        skip_rows : int : default = 0
            number of rows to skip at the top of text files;
            passed to np.loadtxt
            
        average_plddt : bool : False
            average plddt across all Protein objects;
            otherwise uses plddt of the first Protein object
            
        average_bfactor : bool : False
            average bfactor across all Protein objects;
            otherwise uses bfactor of the first Protein object
            
        average_dist_mat : bool : False
            average dist_mat across all Protein objects;
            otherwise uses dist_mat of the first Protein object
            
        """
        # Parameters
        self.max_bfactor =      kwargs.get('max_bfactor', 0.0)
        self.min_plddt =        kwargs.get('min_plddt', 70.0)
        self.neigh_cut =        kwargs.get('neigh_cut', 13.0)
        self.chain =            kwargs.get('chain', '')
        self.fix_pdb =  kwargs.get('fix_pdb', False)
        self.delimiter = kwargs.get('delimiter', ' ')
        self.skip_rows = kwargs.get('skip_rows', 0)

        self.average_plddt =    kwargs.get('average_plddt', False)
        self.average_bfactor =  kwargs.get('average_bfactor', False)
        self.average_dist_mat = kwargs.get('average_bfactor', False)

        # Defaults
        self.plddt = None
        self.bfactor = None
        self.sequence = None
        self.dist_mat = None
        self.neigh_idx = []
        self.neigh_tensor = []
        self.seq_len = None         # Number of residues
        self.num_repeat = None      # Number of structures

        # Load data
        self._parse_input(inp_obj, **kwargs)
        self.get_average_structure()


    def _parse_input(self, inp_obj, **kwargs):
        """
        Parses a list of objects

        Parses a range of input object types (.pdb, .ent, .cif, .txt, .dat, .npy)

        Objects in a list need to be the same type

        Extracts: Protein objects
        """
        self.proteins = []
        # Ensure that input is a list or an array
        if not isinstance(inp_obj, (list, np.ndarray)):
            raise Exception(f"Input object type {type(inp_obj)} is not supported. Acceptable types:" + \
                             "\n\tA list of Protein objects" + \
                             "\n\tA list of Protein paths to pdb / mmcif files" + \
                             "\n\tA list of numpy arrays of atomic coordinates")

        # Parse inputs to get a list of Protein objects
        for item in inp_obj:
            if isinstance(item, Protein):
                self.proteins.append(item)

            elif isinstance(item, (str, Path, np.ndarray)):
                self.proteins.append(Protein(item, **kwargs))

            else:
                raise Exception(f"Input object type {type(item)} is not supported.")

        # Set pLDDT
        if self.average_plddt:
            try:
                self.plddt = np.mean([protein.plddt for protein in self.proteins], axis=0)
            except TypeError:
                self.plddt = None
        else:
            self.plddt = self.proteins[0].plddt

        # Set Bfactor
        if self.average_bfactor:
            try:
                self.bfactor = np.mean([protein.bfactor for protein in self.proteins], axis=0)
            except TypeError:
                self.bfactor = None
        else:
            self.bfactor = self.proteins[0].bfactor

        # Set Distance Matrix
        if self.average_dist_mat:
            self.dist_mat = np.mean([protein.dist_mat for protein in self.proteins], axis=0)
        else:
            self.dist_mat = self.proteins[0].dist_mat.copy()


        self.sequence = self.proteins[0].sequence
        self.seq_len = self.proteins[0].seq_len
        self.num_repeat = len(self.proteins)


    def get_average_structure(self):
        """Get averaged local neighborhoods"""
        self._get_local_neighborhood()

        # Only include neighbors for each residue that are neighbors in all structures
        self._consolidate_neighbor_lists()
        
        # Rotate neighborhood tensors to arbitrary reference neighborhoods,
        # then average over repeat structures
        self._rotate_and_average_neighbor_tensors()


    def _get_local_neighborhood(self):
        """Get local neighborhoods per residue"""
        for protein in self.proteins:
            # If not yet calculated, get the local neighborhood
            if len(self.neigh_idx) == 0:
                kwargs = {a:getattr(self, a) for a in ['min_plddt', 'max_bfactor', 'neigh_cut']}
                for k, v in kwargs.items():
                    setattr(protein, k, v)
                protein.get_local_neighborhood()


    def _consolidate_neighbor_lists(self):
        """Exclude neighbors that are not present in all Protein objects"""
        self.neigh_idx = []
        for i in range(self.seq_len):
            # Count how many times indices are included in neighbor lists
            idx_count = Counter(j for protein in self.proteins for j in protein.neigh_idx[i])
            # Only include indices that are in all neighbor lists
            self.neigh_idx.append([j for j, count in idx_count.items() if count == self.num_repeat])


    ### For each residue j, rotate all neighbourhoods to match the first one,
    ### then average over repeat structures
    def _rotate_and_average_neighbor_tensors(self):
        """Align neighborhood tensors by rotating to match the
           neighborhood tensor of the first Protein object."""
        self.neigh_tensor = []
        for i in range(self.seq_len):
            # If no neighbors..
            if not len(self.neigh_idx[i]):
                self.neigh_tensor.append(np.empty((0,3)))
                continue

            for j in range(self.num_repeat):
                # Only include the rows of the tensor that correspond to
                # the consolidated neighbor list
                idx = [i0 for i0, j0 in enumerate(self.proteins[j].neigh_idx[i]) if j0 in self.neigh_idx[i]]

                if not j:
                    # Do not rotate the first example for residue j
                    # Initialize the list
                    tensor = [self.proteins[j].neigh_tensor[i][idx]]

                else:
                    # Rotate tensor so that it matches the first (reference) tensor
                    rotated_tensor = rotate_points(self.proteins[j].neigh_tensor[i][idx], tensor[0])
                    tensor.append(rotated_tensor)
            self.neigh_tensor.append(np.array(tensor).mean(axis=0))


    def recalculate_average_structure(self):
        """Recalculate average structure after changing parameters""" 
        for protein in self.proteins:
            changed = np.zeros(3, bool)
            for i, attr in enumerate(['min_plddt', 'max_bfactor', 'neigh_cut']):
                old, new = getattr(protein, attr), getattr(self, attr)
                if old != new:
                    setattr(protein, attr, new)
                    changed[i] = True

            # If nothing has changed, exit function
            if np.all(changed == False):
                return

            # If "min_plddt" or "max_bfactor" have changed,
            # recalculate the coords as NaN values may be different
            if np.any(changed[:2] == True):
                protein._update_nan_coords()
                protein._get_dist_mat()

            # Recalculate neighborhoods with updated neighbor cutoff
            protein.get_local_neighborhood()

        self._consolidate_neighbor_lists()
        self._rotate_and_average_neighbor_tensors()

    

