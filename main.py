from pathlib import Path

import numpy as np

from pdb_parser import parse_pdb_coordinates, load_and_fix_pdb_data
from utils.import rotate_points


class Protein:
    def __init__(self, inp_obj, **kwargs):
        # Parameters
        self.max_bfactor = kwargs.get('max_bfactor', 0.0)
        self.min_plddt = kwargs.get('min_plddt', 70.0)
        self.neigh_cut = kwargs.get('neigh_cut', 13.0)
        self.chain = kwargs.get('chain', '')
        self.pdb_fill_missing_nan =  kwargs.get('pdb_fill_missing_nan', False)

        # Defaults
        self.sequence = None
        self.plddt = None
        self.bfactor = None
        self.dist_mat = None
        self.name = kwargs.get('name', '') 
        self.neigh_idx = []
        self.neigh_tensor = []

        # Load data
        self.parse_input(inp_obj)


    def parse_input(self, inp_obj):
        if isinstance(inp_obj, (str, Path)):
            self.load_data_from_path(inp_obj)

        elif isinstance(inp_obj, np.ndarray):
            self.coord = inp_obj
            self.idx = np.arange(len(inp_obj))

        else:
            raise Exception(f"Input object type {type(inp_obj)} is not supported.")


    def load_data_from_path(self, path):
        if self.pdb_fill_missing_nan:
            data = load_and_fix_pdb_data(path, self.chain)
            # (indices in PDB files are inconsistent, so we reorder
            #  the indices so that they start from zero)

        else:
            data = parse_pdb_coordinates(path, self.chain)
            # Reorder indices so that they start from zero
            # (by default, AF-predicted structures start counting at one)
            data[1] = data[1] - 1

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
            self.coord[(self.bfactor > self.max_bfactor)|(np.isnan(bfactor)] = np.nan
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
    def __init__(self, proteins, **kwargs):
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
                self.proteins.append(inp_obj)

            elif isinstance(item, (str, Path, np.ndarray)):
                self.proteins.append(Protein(inp_obj, **kwargs))

            else:
                raise Exception(f"Input object type {type(item)} is not supported.")

        # Set pLDDT
        if self.average_plddt:
            try:
                self.plddt = np.mean([protein.plddt for protein in self.proteins], axis=0)
            except TypeError:
                self.plddt = None
        else:
            self.plddt = proteins[0].plddt

        # Set Bfactor
        if self.average_bfactor:
            try:
                self.bfactor = np.mean([protein.bfactor for protein in self.proteins], axis=0)
            except TypeError:
                self.bfactor = None
        else:
            self.bfactor = proteins[0].bfactor

        # Set Distance Matrix
        if self.average_dist_mat:
            self.dist_mat = np.mean([protein.dist_mat for protein in self.proteins], axis=0)
        else:
            self.dist_mat = proteins[0].dist_mat.copy()


        # Set index
        self.idx = proteins[0].idx
        #    Check that indices are the same, and print a warning if not.
        if self.check_idx_equiv:
            for protein in proteins[1:]:
                if not np.all(self.idx == protein.idx):
                    print("WARNING! Indices are not the same in different proteins." + \
                          "\n\tAre you sure the input list is correct?")

        self.sequence = proteins[0].sequence
        self.seq_len = len(self.idx)
        self.num_repeat = len(self.proteins)


    def get_average_structure(self):
        self._get_local_neighborhood()

        # Only include neighbors for each residue that are neighbors in all structures
        self._consolidate_neighbor_lists()
        
        # Rotate neighborhood tensors to arbitrary reference neighborhoods,
        # then average over repeat structures
        _rotate_and_average_neighbor_tensors()


    def _get_local_neighborhood(self):
        for protein in self.proteins
            # If not yet calculated, get the local neighborhood
            if len(self.neigh_idx) == 0:
                kwargs = {a:getattr(self, a) for a in ['min_plddt', 'max_bfactor', 'neigh_cut']}
                protein.get_local_neighborhood(**kwargs)


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
            for j in range(self.num_repeat):
                # Only include the rows of the tensor that correspond to
                # the consolidated neighbor list
                idx = [i0 for i0, j0 in enumerate(self.proteins[j].neigh_idx[i]) if j0 in self.neigh_idx[j]]

                if not i:
                    # Do not rotate the first example for residue j
                    # Initialize the list
                    tensor = [self.proteins[j].neigh_tensor[i][idx]]

                else:
                    # Rotate tensor so that it matches the first (reference) tensor
                    rotated_tensor = rotate_points(self.proteins[j].neigh_tensor[i][idx], tensor[0])
                    tensor.append(rotated_tensor)
            self.neigh_tensor.append(np.array(tensor).mean(axis=0))


class Deformation:
    def __init__(self, protein_1, protein_2, **kwargs):
        self.prot1 = protein_1
        self.prot2 = protein_2
        self.proteins = [self.prot1, self.prot2]
        
        self.default_method = ["strain"]
        self.all_methods = ["mut_dist", "strain", "shear", "non-affine", "ldd", "lddt", "neighborhood_dist"]
        self.lddt_cutoffs = kwargs.get("lddt_cutoffs", [0.5, 1, 2, 4])
        self.method = kwargs.get('method', self.default_method.copy())
        self.neigh_cut = kwargs.get('neigh_cut', 13.0)
        self.force_cutoff =  kwargs.get('force_cutoff', False) # Recalculate
        self.force_norm =  kwargs.get('force_norm', False)
        self.force_nonorm =  kwargs.get('force_nonorm', False)
        self.force_relative =  kwargs.get('force_relative', False)
        self.force_absolute =  kwargs.get('force_absolute', False)
        self.force_l1norm =  kwargs.get('force_l1norm', False)
        self.force_l2norm =  kwargs.get('force_l2norm', False)
        self.deformation = {}
        self.mutations = []
        self.mutation_idx = []

        self.parse_input()
        self.parse_method()


    def parse_input(self):
        # Check input types
        for prot in self.proteins:
            if not isinstance(self.prot, (AverageProtein, Protein)):
                raise Exception(f"Input object type {type(self.prot)} is not supported.")
            
        # Check that neighbor cutoff definitions are consistent
        _check_neighborhoods()

        # Make sure neighborhoods have been calculated with the correct neigh_cut value
        for prot in self.proteins:
            if isinstance(prot, Protein):


    def _update_protein_neighborhood(self, neigh_cut):
        # If not calculated yet... calculate 
        if len(prot.neigh_idx) == 0:
            prot.neigh_cut = neigh_cut
            prot.get_local_neighborhood()

        # If neigh cut is wrong... calculate 
        print(f"Recalculating Protein with new neighbor cutoff = {neigh_cut}")
        elif prot.neigh_cut != neigh_cut:
            prot.neigh_cut = neigh_cut
            prot.get_local_neighborhood()


    def _update_averageProtein_neighborhood(self, neigh_cut):
        # If not calculated yet... calculate 
        if len(prot.neigh_idx) == 0:
            prot.neigh_cut = neigh_cut
            prot.get_average_structure()

        # If neigh cut is wrong... calculate 
        print(f"Recalculating AverageProtein with new neighbor cutoff = {neigh_cut}")
        elif prot.neigh_cut != neigh_cut:
            prot.neigh_cut = neigh_cut
            prot.recalculate_average_structure()


    # Check for consistency in the use of neighbor cutoffs
    def _check_neighborhoods(self):
        neigh_cut = [prot.neigh_cut for prot in self.proteins]

        # Check for consistency between Protein objects.
        # If they are inconsistent, either force them to use Deformation.neigh_cutoff,
        # or exit
        if neigh_cut[0] != neigh_cut[1]:
            if self.force_cutoff:
                for prot in self.proteins:
                    if isinstance(prot, Protein):
                        self._update_protein_neighborhood(self.neigh_cut)
                    if isinstance(prot, AverageProtein):
                        self._update_averageProtein_neighborhood(self.neigh_cut)
            else:
                raise Exception("AverageProtein / Protein objects were created with different neighbor cutoffs!" + \
                                "\n\tYou need to use the same neighbor cutoff for each structure," + \
                                "\n\tor to automatically recalculate neighborhoods using Deformation.neigh_cutoff," +\
                                " use Deformation(..., force_cutoff=True)")
            return

        # If not "force_cutoff", then ensure Deformation.neigh_cutoff equals that of the Protein objects
        if not self.force_cutoff:
            if self.neigh_cut != neigh_cut[0]:
                self.neigh_cut = neigh_cut[0]
                print(f"WARNING! Resetting neighbour cutoff to {self.neigh_cut}, since this" + \
                       "value was used for the AverageProtein structure." + \
                       "\n\tTo override this, use Deformation(..., force_cutoff=True)")

        for i, prot in enumerate(self.proteins):
            if isinstance(prot, Protein):
                self._update_protein_neighborhood(self.neigh_cut)
            if isinstance(prot, AverageProtein):
                self._update_averageProtein_neighborhood(self.neigh_cut)
            if self.force_cutoff:
                for i, prot in enumerate(self.proteins):
                    prot.neigh_cut = self.neigh_cut
                    prot.recalculate_average_structure()
            else:
                self.neigh_cut = neigh_cut[0]

    
    # Parse method, and ensure methods are acceptable
    def parse_method(self):
        if isinstance(self.method, str):
            if self.method = 'all':
                self.method = self.all_methods.copy()

        elif isinstance(self.method, (list, np.ndarray, set, tuple)):
            method_list = []
            for method in self.method:
                if method in self.all_methods:
                    method_list.append(method)
                else:
                    print(f"WARNING! {method} is not an accepted method")
            self.method = method_list
        else:
            self.method = []

        if len(self.method) == 0:
            raise Exception("No acceptable method found!")
            

    # Set method, and run through parse to check method validity
    def set_method(self, value):
        self.method = value
        self.parse_method()

            
    def calculate_deformation(self):
        # Calculate deformation based on the specified method
        self.deformation = {}
        for method in self.method:
            self._run_analysis(method)


    def _run_analysis(self, method):
        analyses = {'mut_dist': self.calculate_dist_from_mutation(),
                    'strain': self.calculate_strain(),
                    'shear': self.calculate_shear(),
                    'non-affine': self.calculate_non_affine(),
                    'ldd': self.calculate_ldd(),
                    'lddt': self.calculate_lddt(),
                    'neighborhood_dist': self.calculate_neighborhood_dist()}


    # Calculate distance from closest mutation
    def calculate_dist_from_mutation(self):
        try:
            self.sub_pos = utils.get_mutation_position(self.prot1.sequence, self.prot2.sequence)
        except AttributeError as E:
            raise AttributeError("Sequence is not defined for Protein object")

        # If none differ, then return np.nan
        if not len(self.sub_pos):
            print("WARNING! Trying to calculate distance from mutation, while comparing identical sequences")
            return np.zeros(self.prot1.seq_len) * np.nan
        
        # Calculate mindist using the full array (inc. nan)
        mut_dist1 = self.prot1.dist_mat[:,sub_pos]
        mut_dist2 = self.prot2.dist_mat[:,sub_pos]

        # Average the distance across both structures,
        # and get the minimum distance per residue to a mutated position
        self.mut_dist = np.nanmin(0.5 * (mut_dist1 + mut_dist2), axis=1)


    def calculate_lddt():
        lddt = np.zeros(self.prot1.seq_len, float) * np.nan
        for i in range(self.prot1.seq_len):
            # If no data for residue, leave np.nan
            if (not len(self.prot1.neigh_idx[i])) | (not len(self.prot2.neigh_idx[i])):
                continue

            # Get shared indices
            i1, i2 = utils.get_shared_indices(self.prot1.neigh_idx[i], self.prot2.neigh_idx[i])

            # If no shared indices, leave np.nan
            if not len(i1):
                continue

            # Get local distance vectors
            v1 = np.linalg.norm(self.prot1.neigh_tensor[i][i1], axis=1)
            v2 = np.linalg.norm(self.prot2.neigh_tensor[i][i2], axis=1)

            # Get local distance difference vector
            dv = v2 - v1

            lddt[i] = np.sum([np.sum(dv<=cut) for cut in self.lddt_cutoffs]) / (len(lddt_cutoffs) * len(dv))

        self.lddt = lddt


    def calculate_ldd(self):
        ldd = np.zeros(self.prot1.seq_len, float) * np.nan
        for i in range(self.prot1.seq_len):
            # If no data for residue, leave np.nan
            if (not len(self.prot1.neigh_idx[i])) | (not len(self.prot2.neigh_idx[i])):
                continue

            # Get shared indices
            i1, i2 = utils.get_shared_indices(self.prot1.neigh_idx[i], self.prot2.neigh_idx[i])

            # If no shared indices, leave np.nan
            if not len(i1):
                continue

            # Get local distance vectors
            v1 = np.linalg.norm(self.prot1.neigh_tensor[i][i1], axis=1)
            v2 = np.linalg.norm(self.prot2.neigh_tensor[i][i2], axis=1)

            # Get local distance difference vector
            dv = v2 - v1

            if self.force_relative:
                dv = dv / v1

            # Calculate local distance difference 
            if self.force_norm:
                ldd[i] = np.linalg.norm(ld1 - ld2) / len(i1)
            else:
                ldd[i] = np.linalg.norm(ld1 - ld2)

        self.ldd = ldd


    def calculate_neighborhood_distance():
        nd = np.zeros(self.prot1.seq_len, float) * np.nan
        for i in range(self.prot1.seq_len):
            # If no data for residue, leave np.nan
            if (not len(self.prot1.neigh_idx[i])) | (not len(self.prot2.neigh_idx[i])):
                continue

            # Get shared indices
            i1, i2 = utils.get_shared_indices(self.prot1.neigh_idx[i], self.prot2.neigh_idx[i])

            # If no shared indices, leave np.nan
            if not len(i1):
                continue

            # Get neighbourhood tensors
            c1 = self.prot1.neigh_tensor[i][i1]
            c2 = self.prot2.neigh_tensor[i][i2]

            # Rotate neighbourhood tensor and calculate Euclidean distance
            if self.force_norm:
                nd[i] = np.linalg.norm(AS.rotate_points(c2, c1) - c1) / len(i1)
            else:
                nd[i] = np.linalg.norm(AS.rotate_points(c2, c1) - c1)

        self.neighbor_distance = nd


    def calculate_shear():
        shear = np.zeros(self.prot1.seq_len, float) * np.nan
        for i in range(self.prot1.seq_len):
            # If no data for residue, leave np.nan
            if (not len(self.prot1.neigh_idx[i])) | (not len(self.prot2.neigh_idx[i])):
                continue

            # Get shared indices
            i1, i2 = utils.get_shared_indices(self.prot1.neigh_idx[i], self.prot2.neigh_idx[i])

            # If no shared indices, leave np.nan
            if not len(i1):
                continue

            # Get neighbourhood tensors
            u1 = self.prot1.neigh_tensor[i][i1]
            u2 = self.prot2.neigh_tensor[i][i2]
            try:
                duu = u2 @ u2.T - u1 @ u1.T
                uu = np.linalg.inv(u1.T @ u1) 
                C = 0.5 * (uu @ u1.T @ duu @ u1 @ uu)
                shear[i] 0.5 * np.sum(np.diag(C@C) - np.diag(C)**2)
            except Exception as e:
                continue

        self.shear = shear


    def calculate_strain():
        strain = np.zeros(self.prot1.seq_len, float) * np.nan
        for i in range(self.prot1.seq_len):
            # If no data for residue, leave np.nan
            if (not len(self.prot1.neigh_idx[i])) | (not len(self.prot2.neigh_idx[i])):
                strain.append(np.nan)
                continue

            # Get shared indices
            i1 = [j for j, k in enumerate(self.prot1.neigh_idx[i]) if k in self.prot2.neigh_idx[i]]
            i2 = [j for j, k in enumerate(self.prot2.neigh_idx[i]) if k in self.prot1.neigh_idx[i]]
            # If no shared indices, leave np.nan
            if not len(i1):
                strain.append(np.nan)
                continue

            # Get neighbourhood tensors
            c1 = self.prot1.neigh_tensor[i][i1]
            c2 = self.prot2.neigh_tensor[i][i2]

            # Rotate neighbourhood tensor and calculate Euclidean distance
            c3 = AS.rotate_points(c2, c1)
            s = 0.0
            for j in range(len(c1)):
                if self.force_absolute:
                    s += np.linalg.norm(c3[j] - c1[j])
                else:
                    s += np.linalg.norm(c3[j] - c1[j]) / np.linalg.norm(c1[j])

            if not self.force_nonorm:
                s /= len(i1)

            strain[i] = s

        return np.array(strain)

