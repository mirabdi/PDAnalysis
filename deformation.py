from pathlib import Path

from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np
import pandas as pd

from protein import Protein, AverageProtein
from utils import rotate_points, get_shared_indices, get_mutation_position


class Deformation:
    def __init__(self, protein_1, protein_2, **kwargs):
        self.prot1 = protein_1
        self.prot2 = protein_2
        self.proteins = [self.prot1, self.prot2]
        
        self.default_method = ["strain"]
        self.all_methods = ["mut_dist", "strain", "shear", "non-affine", "ldd", "lddt", "neighborhood_dist", "rmsd"]
        self.lddt_cutoffs = kwargs.get("lddt_cutoffs", [0.5, 1, 2, 4])
        self.method = kwargs.get('method', self.default_method.copy())
        self.neigh_cut = kwargs.get('neigh_cut', 13.0)

        self.force_cutoff =  kwargs.get('force_cutoff', False)
        self.force_norm =  kwargs.get('force_norm', False)
        self.force_nonorm =  kwargs.get('force_nonorm', False)
        self.force_relative =  kwargs.get('force_relative', False)
        self.force_absolute =  kwargs.get('force_absolute', False)
        self.force_l1norm =  kwargs.get('force_l1norm', False)
        self.force_l2norm =  kwargs.get('force_l2norm', False)

        self.verbose = kwargs.get('verbose', False)

        self.deformation = {}
        self.mutations = []
        self.mutation_idx = []

        self.parse_input()
        self.parse_method()

        if self.verbose:
            self.print_inputs_summary()


    def print_inputs_summary(self):
        print(f"Comparing {self.prot1} with {self.prot2}.")
        print(f"Sequence length :: {self.prot1.seq_len}")

        nmiss1 = sum([len(i) == 0 for i in self.prot1.neigh_idx])
        nmiss2 = sum([len(i) == 0 for i in self.prot2.neigh_idx])
        print(f"Number of residues excluded due to missing coordinates, or due to low pLDDT" + \
              f" / high B-factor ::\n\tProtA, {nmiss1}\n\tProtB, {nmiss2}")

        print(f"Amino acid substitutions :: {' '.join(self.sub_str)}")
        print(f"Methods to run :: {' '.join(self.method)}")


    def parse_input(self):
        # Check input types
        for prot in self.proteins:
            if not isinstance(prot, (AverageProtein, Protein)):
                raise Exception(f"Input object type {type(prot)} is not supported.")
            
        # Check that neighbor cutoff definitions are consistent.
        # This method also loads neighborhoods if they are not loaded already.
        self._check_neighborhoods()

        # Check that protein coordinate arrays are the same length
        l1 = len(self.prot1.neigh_idx)
        l2 = len(self.prot2.neigh_idx)
        if l1 != l2:
            raise Exception("Protein coordinate arrays are not the same length: " + \
                  f"Protein A has {l1} residues, while " + \
                  f"Protein B has {l2} residues.\n" + \
                  "If using PDB files with missing coordinates, use the --pdb_fill_missing_nan option.")

        try:
            self.sub_pos = get_mutation_position(self.prot1.sequence, self.prot2.sequence)
            self.sub_str = self.get_substitution_strings()
        except AttributeError as E:
            raise AttributeError("Sequence is not defined for Protein object")


    def get_substitution_strings(self):
        return [f"{self.prot1.sequence[i]}{i+1}{self.prot2.sequence[i]}" for i in self.sub_pos]


    def _update_protein_neighborhood(self, prot, neigh_cut):
        # If not calculated yet... calculate 
        if len(prot.neigh_idx) == 0:
            prot.neigh_cut = neigh_cut
            prot.get_local_neighborhood()

        # If neigh cut is wrong... calculate 
        elif prot.neigh_cut != neigh_cut:
            print(f"Recalculating Protein with new neighbor cutoff = {neigh_cut}")
            prot.neigh_cut = neigh_cut
            prot.get_local_neighborhood()


    def _update_averageProtein_neighborhood(self, prot, neigh_cut):
        # If not calculated yet... calculate 
        if len(prot.neigh_idx) == 0:
            prot.neigh_cut = neigh_cut
            prot.get_average_structure()

        # If neigh cut is wrong... calculate 
        elif prot.neigh_cut != neigh_cut:
            print(f"Recalculating AverageProtein with new neighbor cutoff = {neigh_cut}")
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
                        self._update_protein_neighborhood(prot, self.neigh_cut)
                    if isinstance(prot, AverageProtein):
                        self._update_averageProtein_neighborhood(prot, self.neigh_cut)
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
                print(f"WARNING! Resetting neighbour cutoff to {self.neigh_cut}, since this " + \
                       "value was used for the AverageProtein structure." + \
                       "\n\tTo override this, use Deformation(..., force_cutoff=True)")

        for i, prot in enumerate(self.proteins):
            if isinstance(prot, Protein):
                self._update_protein_neighborhood(prot, self.neigh_cut)
            if isinstance(prot, AverageProtein):
                self._update_averageProtein_neighborhood(prot, self.neigh_cut)
            if self.force_cutoff:
                for i, prot in enumerate(self.proteins):
                    prot.neigh_cut = self.neigh_cut
                    prot.recalculate_average_structure()
            else:
                self.neigh_cut = neigh_cut[0]

    
    # Parse method, and ensure methods are acceptable
    def parse_method(self):
        if isinstance(self.method, (list, np.ndarray, set, tuple)):
            if len(self.method) == 1:
                if self.method[0] == 'all':
                    self.method = self.all_methods.copy()

            else:
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
            raise Exception("No acceptable method found!" + \
            f"\nChoose one out of: all, {', '.join(self.all_methods)}")
            

    # Set method, and run through parse to check method validity
    def set_method(self, value):
        self.method = value
        self.parse_method()


    def save_output(self, path_out):
        # Load any deformation that was calculated
        deform = {}
        for m in self.all_methods:
            if hasattr(self, m):
                if m != 'rmsd':
                    deform[m] = getattr(self, m)
                else:
                    # Since RMSD is a scalar, we create a list to match the output format
                    deform[m] = [getattr(self, m)] * self.prot1.seq_len

        # Only save output if deformation was calculated
        if not len(deform):
            raise Exception("ERROR! Cannot save output if there is none!")

        # Add sequence indices, residue names
        else:
            type_names = ['residue_index', 'protA_resname', 'protB_resname']
            res_data = [np.arange(len(self.prot1.sequence)) + 1, self.prot1.sequence, self.prot2.sequence]
            output = {k: d for k, d in zip(type_names, res_data) if not isinstance(d, type(None))}
            ########################################
            ### NEED TO ADD NUMBER OF NEIGHBORS
            ########################################

            output.update(deform)
    
            pd.DataFrame(output).to_csv(path_out, index=False)

            
    def calculate_deformation(self):
        # Calculate deformation based on the specified method
        self.deformation = {}
        for method in self.method:
            self._run_analysis(method)


    def _run_analysis(self, method):
        analyses = {'mut_dist': self.calculate_dist_from_mutation,
                    'strain': self.calculate_strain,
                    'shear': self.calculate_shear,
                    'non-affine': self.calculate_non_affine,
                    'ldd': self.calculate_ldd,
                    'lddt': self.calculate_lddt,
                    'neighborhood_dist': self.calculate_neighborhood_dist,
                    'rmsd': self.calculate_rmsd}
        analyses[method]()


    def _get_shared_indices(self):
        self.shared_indices = []
        for i in range(self.prot1.seq_len):
            # If no data for residue, shared indices is empty
            if (not len(self.prot1.neigh_idx[i])) | (not len(self.prot2.neigh_idx[i])):
                self.shared_indices.append(([], []))
            else:
                self.shared_indices.append(get_shared_indices(self.prot1.neigh_idx[i], self.prot2.neigh_idx[i]))


    # Calculate distance from closest mutation
    def calculate_dist_from_mutation(self):
        # If none differ, then return np.nan
        if not len(self.sub_pos):
            print("WARNING! Trying to calculate distance from mutation, while comparing identical sequences")
            return np.zeros(self.prot1.seq_len) * np.nan
        
        # Calculate mindist using the full array (inc. nan)
        mut_dist1 = self.prot1.dist_mat[:,self.sub_pos]
        mut_dist2 = self.prot2.dist_mat[:,self.sub_pos]

        # Average the distance across both structures,
        # and get the minimum distance per residue to a mutated position
        self.mut_dist = np.nanmin(0.5 * (mut_dist1 + mut_dist2), axis=1)


    def calculate_lddt(self):
        lddt = np.zeros(self.prot1.seq_len, float) * np.nan

        if not hasattr(self, "shared_indices"):
            self._get_shared_indices()

        for i in range(self.prot1.seq_len):
            # Get shared indices
            i1, i2 = self.shared_indices[i]

            # If no shared indices, leave np.nan
            if not len(i1):
                continue

            # Get local distance vectors
            v1 = np.linalg.norm(self.prot1.neigh_tensor[i][i1], axis=1)
            v2 = np.linalg.norm(self.prot2.neigh_tensor[i][i2], axis=1)

            # Get local distance difference vector
            dv = v2 - v1

            lddt[i] = np.sum([np.sum(dv<=cut) for cut in self.lddt_cutoffs]) / (len(self.lddt_cutoffs) * len(dv))

        self.lddt = lddt


    def calculate_ldd(self):
        ldd = np.zeros(self.prot1.seq_len, float) * np.nan

        if not hasattr(self, "shared_indices"):
            self._get_shared_indices()

        for i in range(self.prot1.seq_len):
            # Get shared indices
            i1, i2 = self.shared_indices[i]

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
                ldd[i] = np.linalg.norm(dv) / len(i1)
            else:
                ldd[i] = np.linalg.norm(dv)

        self.ldd = ldd


    def calculate_neighborhood_dist(self):
        nd = np.zeros(self.prot1.seq_len, float) * np.nan

        if not hasattr(self, "shared_indices"):
            self._get_shared_indices()

        for i in range(self.prot1.seq_len):
            # Get shared indices
            i1, i2 = self.shared_indices[i]

            # If no shared indices, leave np.nan
            if not len(i1):
                continue

            # Get neighbourhood tensors
            c1 = self.prot1.neigh_tensor[i][i1]
            c2 = self.prot2.neigh_tensor[i][i2]

            # Rotate neighbourhood tensor and calculate Euclidean distance
            if self.force_norm:
                nd[i] = np.linalg.norm(rotate_points(c2, c1) - c1) / len(i1)
            else:
                nd[i] = np.linalg.norm(rotate_points(c2, c1) - c1)

        self.neighbor_distance = nd


    def calculate_shear(self):
        shear = np.zeros(self.prot1.seq_len, float) * np.nan

        if not hasattr(self, "shared_indices"):
            self._get_shared_indices()

        for i in range(self.prot1.seq_len):
            # Get shared indices
            i1, i2 = self.shared_indices[i]

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
                shear[i] = 0.5 * np.sum(np.diag(C@C) - np.diag(C)**2)
            except Exception as e:
                continue

        self.shear = shear


    def calculate_strain(self):
        strain = np.zeros(self.prot1.seq_len, float) * np.nan

        if not hasattr(self, "shared_indices"):
            self._get_shared_indices()

        for i in range(self.prot1.seq_len):
            # Get shared indices
            i1, i2 = self.shared_indices[i]

            # If no shared indices, leave np.nan
            if not len(i1):
                continue

            # Get neighbourhood tensors
            c1 = self.prot1.neigh_tensor[i][i1]
            c2 = self.prot2.neigh_tensor[i][i2]

            # Rotate neighbourhood tensor and calculate Euclidean distance
            c3 = rotate_points(c2, c1)
            if self.force_absolute:
                s = np.sum(np.linalg.norm(c3 - c1, axis=1))
            else:
                s = np.sum(np.linalg.norm(c3 - c1, axis=1) / np.linalg.norm(c1, axis=1))

            if not self.force_nonorm:
                s /= len(i1)

            strain[i] = s

        self.strain = strain

    def calculate_non_affine(self):
        pass

    #######################################################
    ### Root-Mean-Square Deviation
    ### Superimpose structures and calculate RMSD
    def calculate_rmsd(self):
        # Cannot calculate RMSD with averaged neighborhoods,
        # so if an AverageProtein object is passsed, we simply
        # take the first Protein object 
        coords = []
        for prot in self.proteins:
            if isinstance(prot, Protein):
                coords.append(prot.coord)
            elif isinstance(prot, AverageProtein):
                coords.append(prot.proteins[0].coord)
                print("WARNING! Trying to calculate RMSD with an AverageProtein object. " + \
                      "Since this is not possible, an embedded Protein object is used instead.")
        c1, c2 = coords

        # Remove NAN values
        idx = ~np.any(np.isnan(c1) | np.isnan(c2), axis=1)
        c1, c2 = c1[idx], c2[idx]

        sup = SVDSuperimposer()
        sup.set(c1, c2) 
        sup.run()
        self.rmsd = sup.get_rms()


