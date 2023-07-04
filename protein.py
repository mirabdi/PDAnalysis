from collections import Counter
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

from pdb_parser import parse_pdb_coordinates, parse_mmcif_coordinates, load_and_fix_pdb_data
from utils import rotate_points


class Protein:
    def __init__(self, inp_obj, **kwargs):
        # Parameters
        self.max_bfactor = kwargs.get('max_bfactor', 0.0)
        self.min_plddt = kwargs.get('min_plddt', 0.0)
        self.neigh_cut = kwargs.get('neigh_cut', 13.0)
        self.chain = kwargs.get('chain', '')
        self.pdb_fill_missing_nan =  kwargs.get('pdb_fill_missing_nan', False)

        # Defaults
        self.sequence = None
        self.idx = None
        self.plddt = None
        self.bfactor = None
        self.dist_mat = None
        self.name = kwargs.get('name', '') 
        self.file_fmt = kwargs.get('file_fmt', '') 
        self.delimiter = kwargs.get('delimiter', ' ')
        self.skip_rows = kwargs.get('skip_rows', 0)
        self.neigh_idx = []
        self.neigh_tensor = []
        self.seq_len = None         # Number of residues

        # Load data
        self.parse_input(inp_obj)
        self.get_local_neighborhood()


    def parse_input(self, inp_obj):
        if isinstance(inp_obj, (str, Path)):
            ext = Path(inp_obj).suffix

            ### Read PDB / mmCIF file
            if (ext in ['.pdb', '.cif']) or (self.file_fmt == ['pdb', 'mmcif', 'cif']):
                self.load_data_from_path(inp_obj, ext)

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

        self.seq_len = len(self.idx)


    def load_data_from_path(self, path, ext):
        if self.pdb_fill_missing_nan:
            data = load_and_fix_pdb_data(path, self.chain)
            # (indices in PDB files are inconsistent, so we reorder
            #  the indices so that they start from zero)

        else:
            if ext == ".pdb":
                data = parse_pdb_coordinates(path, self.chain)
            elif ext == ".cif":
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
        self.coord = self.coord_raw.copy()
        if self.max_bfactor > 0.0:
            self.coord[(self.bfactor > self.max_bfactor)|(np.isnan(bfactor))] = np.nan
        elif self.min_plddt > 0.0:
            self.coord[self.plddt <= self.min_plddt] = np.nan


    def _get_dist_mat(self):
        self.dist_mat = cdist(self.coord, self.coord)
    

    def get_local_neighborhood(self):
        if not isinstance(self.dist_mat, np.ndarray):
            self._get_dist_mat()
        self.neigh_idx = [np.where((d > 0) & (d <= self.neigh_cut) & (np.isfinite(d)))[0] for d in self.dist_mat]
        self._calculate_neighbor_tensor()


    def _calculate_neighbor_tensor(self):
        self.neigh_tensor = []
        for i, idx in enumerate(self.neigh_idx):
            self.neigh_tensor.append(self.coord[idx] - self.coord[[i]])
        
    



class AverageProtein:
    def __init__(self, inp_obj, **kwargs):
        # Parameters
        self.max_bfactor =      kwargs.get('max_bfactor', 0.0)
        self.min_plddt =        kwargs.get('min_plddt', 70.0)
        self.neigh_cut =        kwargs.get('neigh_cut', 13.0)
        self.chain =            kwargs.get('chain', '')
        self.average_plddt =    kwargs.get('average_plddt', False)
        self.average_bfactor =  kwargs.get('average_bfactor', False)
        self.average_dist_mat = kwargs.get('average_bfactor', False)
        self.check_idx_equiv =  kwargs.get('force_idx_equiv', True)
        self.pdb_fill_missing_nan =  kwargs.get('pdb_fill_missing_nan', False)

        # Defaults
        self.plddt = None
        self.bfactor = None
        self.idx = None
        self.sequence = None
        self.dist_mat = None
        self.name = kwargs.get('name', '')
        self.neigh_idx = []
        self.neigh_tensor = []
        self.neigh_tensor = []
        self.seq_len = None         # Number of residues
        self.num_repeat = None      # Number of structures

        # Load data
        self.parse_input(inp_obj, **kwargs)
        self.get_average_structure()


    def parse_input(self, inp_obj, **kwargs):
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


        # Set index
        self.idx = self.proteins[0].idx
        #    Check that indices are the same, and print a warning if not.
        if self.check_idx_equiv:
            for protein in self.proteins[1:]:
                if not np.all(self.idx == protein.idx):
                    print("WARNING! Indices are not the same in different proteins." + \
                          "\n\tAre you sure the input list is correct?")

        self.sequence = self.proteins[0].sequence
        self.seq_len = self.proteins[0].seq_len
        self.num_repeat = len(self.proteins)


    def get_average_structure(self):
        self._get_local_neighborhood()

        # Only include neighbors for each residue that are neighbors in all structures
        self._consolidate_neighbor_lists()
        
        # Rotate neighborhood tensors to arbitrary reference neighborhoods,
        # then average over repeat structures
        self._rotate_and_average_neighbor_tensors()


    def _get_local_neighborhood(self):
        for protein in self.proteins:
            # If not yet calculated, get the local neighborhood
            if len(self.neigh_idx) == 0:
                kwargs = {a:getattr(self, a) for a in ['min_plddt', 'max_bfactor', 'neigh_cut']}
                for k, v in kwargs.items():
                    setattr(protein, k, v)
                protein.get_local_neighborhood()


    def recalculate_average_structure(self):
        for protein in self.proteins:
            changed = np.zeros(3, bool)
            for i, attr in ['min_plddt', 'max_bfactor', 'neigh_cut']:
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

            # Recalculate neighborhoods with updated neighbor cutoff
            protein.get_local_neighborhood()

    
    def _consolidate_neighbor_lists(self):
        self.neigh_idx = []
        for i in range(self.seq_len):
            # Count how many times indices are included in neighbor lists
            idx_count = Counter(j for protein in self.proteins for j in protein.neigh_idx[i])
            # Only include indices that are in all neighbor lists
            self.neigh_idx.append([j for j, count in idx_count.items() if count == self.num_repeat])


    ### For each residue j, rotate all neighbourhoods to match the first one,
    ### then average over repeat structures
    def _rotate_and_average_neighbor_tensors(self):
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


