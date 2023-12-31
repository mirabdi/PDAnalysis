from pathlib import Path

import numpy as np
import pandas as pd

from .protein import Protein, AverageProtein
from .utils import rotate_points, get_shared_indices, get_mutation_position


class Deformation:
    """ Deformation object description:

    The Deformation object takes two protein objects (Protein, AverageProtein),
    and calculates a range of deformation metrics.

    Methods
    -------

    Effective Strain : strain
        mean relative change in distance between neighbors in two neighborhoods
        [McBride, J. M., Polev, K., Abdirasulov, A., Reinharz, V., Grzybowski, B. A.,
         & Tlusty, T. (2023), 'AlphaFold2 can predict single-mutation effects', biorXiv]

    Shear Strain : shear
        magnitude of the off-diagonal components of the strain tensor
        [Eckmann, J P, J. Rougemont, and T. Tlusty (2019), “Colloquium: Proteins:
         The physics of amorphous evolving matter,” Rev. Mod. Phys. 91, 031001]
        
    Non-Affine Strain : non_affine
        non-linear component of strain, defined by the residual of the fit of the strain tensor
        to two neighborhood tensors
        [Falk, M L, and J. S. Langer (1998), 'Dynamics of viscoplastic deformation in
         amorphous solids', Phys. Rev. E 57, 7192–7205]

    Local Distance Difference (LDD) : ldd
        L2 norm of the difference between distances between neighbors in two neighborhoods
        [McBride, J. M., Polev, K., Abdirasulov, A., Reinharz, V., Grzybowski, B. A.,
         & Tlusty, T. (2023), 'AlphaFold2 can predict single-mutation effects', biorXiv]

    Local Distance Difference Test (LDDT) : lddt
        fraction of differences between distances between neighbors that are within some
        set thresholds in two neighborhoods
        [Mariani, Valerio, Marco Biasini, Alessandro Barbato, and Torsten Schwede (2013),
         'lDDT: a local superposition-free score for comparing protein structures and models
         using distance difference tests', Bioinformatics 29 (21), 2722–2728]

    Neighborhood Distance : neighborhood_dist
        L2 norm of the difference between two neighborhood tensors
        [McBride, J. M., Polev, K., Abdirasulov, A., Reinharz, V., Grzybowski, B. A.,
         & Tlusty, T. (2023), 'AlphaFold2 can predict single-mutation effects', biorXiv]

    Root Mean Square Deviation (RMSD) : rmsd
        L2 norm of the difference between two protein structures

    Distance to Mutated Site : mut_dist
        Distance (angstroms) of each residue to the nearest mutated site


    Attributes
    ----------

    all_methods : list
        list of all methods:
            ['mut_dist', 'strain', 'shear', 'non_affine', 'ldd',
             'lddt', 'neighborhood_dist', 'rmsd']

    sub_pos : list(int)
        list of sequence indices where sequences of protein_1 and protein_2 are different;
        indices start from zero

    sub_str : list(str)
        sequence substitutions written in the format {amino_acid_1}{index}{amino_acid_2};
        indices start from one (as is convention)

    """
    def __init__(self, protein_1, protein_2, **kwargs):
        """
        args
        ----------
        protein_1: (Protein, AverageProtein)
        protein_2: (Protein, AverageProtein)


        kwargs
        ----------

        method : (str, list, set, tuple, np.ndarray) : default = ["strain"]
            method(s) to be used when calling self.run();
            'all' results in all methods being used;
            multiple methods can be passed as an iterable object

        lddt_cutoffs : list : default = [0.5, 1, 2, 4]
            specify a list of thresholds used to calculate LDDT

        neigh_cut : float : default = 13.0
            cutoff radius used to determine neighborhood;
            can be used to update Protein and AverageProtein objects

        force_cutoff : bool : default = False
            recalculate Protein and AverageProtein neighborhoods if
            Deformation.neigh_cut != Protein.neigh_cut

        force_norm : bool : default = False
            force a metric that is not normally normalized by the number
            of neighbors to be normalized in this way

        force_nonorm : bool : default = False
            force a metric that is normally normalized by the number
            of neighbors to not be normalized in this way

        force_relative : bool : default = False
            force a metric that is not normally normalized by neighbor
            distance to be normalized in this way

        force_absolute : bool : default = False
            force a metric that is normally normalized by neighbor
            distance to not be normalized in this way

        force_nonorm : bool : default = False
            force a metric that is normally normalized by the number
            of neighbors to not be normalized in this way

        verbose : bool : default = True
            print a summary of the deformation calculation

        """
        self.prot1 = protein_1
        self.prot2 = protein_2
        self.proteins = [self.prot1, self.prot2]
        
        self.default_method = ["strain"]
        self.all_methods = ["mut_dist", "strain", "shear", "non_affine", "ldd", "lddt", "neighborhood_dist", "rmsd"]
        self.lddt_cutoffs = kwargs.get("lddt_cutoffs", [0.5, 1, 2, 4])
        self.method = kwargs.get('method', self.default_method.copy())
        self.neigh_cut = kwargs.get('neigh_cut', 13.0)

        self.force_cutoff =  kwargs.get('force_cutoff', False)
        self.force_norm =  kwargs.get('force_norm', False)
        self.force_nonorm =  kwargs.get('force_nonorm', False)
        self.force_relative =  kwargs.get('force_relative', False)
        self.force_absolute =  kwargs.get('force_absolute', False)

        self.verbose = kwargs.get('verbose', True)

        self.sub_pos = None
        self.sub_str = ''

        self._parse_input()
        self._parse_method()

        if self.verbose:
            self._print_inputs_summary()


    def _print_inputs_summary(self):
        """Print summary of the deformation calculation"""
        print(f"Comparing {self.prot1} with {self.prot2}.")
        print(f"Sequence length :: {self.prot1.seq_len}")

        nmiss1 = sum([len(i) == 0 for i in self.prot1.neigh_idx])
        nmiss2 = sum([len(i) == 0 for i in self.prot2.neigh_idx])
        print(f"Number of residues excluded due to missing coordinates, or due to low pLDDT" + \
              f" / high B-factor ::\n\tProtA, {nmiss1}\n\tProtB, {nmiss2}")

        print(f"Amino acid substitutions :: {' '.join(self.sub_str)}")
        print(f"Methods to run :: {' '.join(self.method)}")


    def _parse_input(self):
        """
        Accepts protein (Protein, AverageProtein) objects as input.

        Checks protein objects for consistency in neighborhood cutoff radii,
        protein size, and checks sequences for differences (mutations).
        """
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
                  "If using PDB files with missing coordinates, use the --fix_pdb option.")

        try:
            self.sub_pos = get_mutation_position(self.prot1.sequence, self.prot2.sequence)
            self.sub_str = self._get_substitution_strings()
        except AttributeError as E:
            raise AttributeError("Sequence is not defined for Protein object")


    def _get_substitution_strings(self):
        """Get conventional representation of mutation as a string"""
        return [f"{self.prot1.sequence[i]}{i+1}{self.prot2.sequence[i]}" for i in self.sub_pos]


    def _update_protein_neighborhood(self, prot, neigh_cut):
        """Check Protein neighbood and recalculate if neigh_cut is wrong"""
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
        """Check AverageProtein neighbood and recalculate if neigh_cut is wrong"""
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
        """
        Run some consistency checks on the two (Protein, AverageProtein) objects.

        Recalculate neighborhoods if parameters (neigh_cut, min_plddt, max_bfactor) have changed.
        """
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
    def _parse_method(self):
        """Parse method into the appropriate format"""
        if isinstance(self.method, str):
            self.method = [self.method]

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
        """Set the method(s) to be used by self.run()
            > (str, list, set, tuple, np.ndarray)
                > 'all' results in all methods being used;
                > multiple methods can be passed as an iterable object
            > list of all methods:
                > 'mut_dist'
                > 'strain'
                > 'shear'
                > 'non_affine'
                > 'ldd'
                > 'lddt'
                > 'neighborhood_dist'
                > 'rmsd'
        """
        self.method = value
        self._parse_method()


    def save_output(self, path_out):
        """Save outputs to CSV file"""
        # Load any deformation that was calculated
        deform = {}
        for m in self.all_methods:
            if hasattr(self, m):
                if m != 'rmsd':
                    deform[m] = getattr(self, m)
                else:
                    deform["rmsd_per_residue"] = self.rmsd_per_residue
                    deform["rmsd_overall"] = [self.rmsd] * len(self.rmsd_per_residue)
                    print(deform["rmsd_per_residue"].shape)

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

            
    def run(self):
        """
        Runs through all of the methods in self.method.

        The output of each method, {m}, is stored in self.{m} (e.g. self.strain);
        RMSD is stored as self.rmsd (whole protein), and as RMSD per residue, self.rmsd_per_residue
        """
        # Calculate deformation based on the specified method
        for method in self.method:
            self._run_analysis(method)


    def _run_analysis(self, method):
        """Run the analysis code for a particular method"""
        analyses = {'mut_dist': self.calculate_dist_from_mutation,
                    'strain': self.calculate_strain,
                    'shear': self.calculate_shear,
                    'non_affine': self.calculate_non_affine,
                    'ldd': self.calculate_ldd,
                    'lddt': self.calculate_lddt,
                    'neighborhood_dist': self.calculate_neighborhood_dist,
                    'rmsd': self.calculate_rmsd}
        analyses[method]()


    def _get_shared_indices(self):
        """Get shared indices between two neighborhoods"""
        self.shared_indices = []
        for i in range(self.prot1.seq_len):
            # If no data for residue, shared indices is empty
            if (not len(self.prot1.neigh_idx[i])) | (not len(self.prot2.neigh_idx[i])):
                self.shared_indices.append(([], []))
            else:
                self.shared_indices.append(get_shared_indices(self.prot1.neigh_idx[i], self.prot2.neigh_idx[i]))


    # Calculate distance from closest mutation
    def calculate_dist_from_mutation(self):
        """Calcualte distance from the nearest mutated residue"""
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


    def _calculate_deformation(self, deformation_method):
        """Find shared residue indices between two neighborhoods
        and calculate deformation per residue"""
        deformation = np.zeros(self.prot1.seq_len, float) * np.nan
        if not hasattr(self, "shared_indices"):
            self._get_shared_indices()

        kwargs = {arg: getattr(self, arg) for arg in ["force_relative", "force_norm", "force_absolute", "force_nonorm"]}
        for i in range(self.prot1.seq_len):
            # Get shared indices
            i1, i2 = self.shared_indices[i]

            # If no shared indices, leave np.nan
            if not len(i1):
                continue

            deformation[i] = deformation_method(self.prot1.neigh_tensor[i][i1], self.prot2.neigh_tensor[i][i2], **kwargs)

        return deformation


    def _calculate_lddt_residue(self, neigh_tensor1, neigh_tensor2, **kwargs):
        """Calculate LDDT given a pair of neighborhood tensors"""
        # Get local distance vectors
        v1 = np.linalg.norm(neigh_tensor1, axis=1)
        v2 = np.linalg.norm(neigh_tensor2, axis=1)

        # Get local distance difference vector
        dv = v2 - v1

        return np.sum([np.sum(dv<=cut) for cut in self.lddt_cutoffs]) / (len(self.lddt_cutoffs) * len(dv))

    def calculate_lddt(self):
        """Calculate LDDT"""
        self.lddt = self._calculate_deformation(self._calculate_lddt_residue)


    def _calculate_ldd_residue(self, neigh_tensor1, neigh_tensor2, **kwargs):
        """Calculate LDD given a pair of neighborhood tensors"""
        # Get local distance vectors
        v1 = np.linalg.norm(neigh_tensor1, axis=1)
        v2 = np.linalg.norm(neigh_tensor2, axis=1)

        # Get local distance difference vector
        dv = v2 - v1

        if kwargs["force_relative"]:
            # Normalize LDD by distance
            dv = dv / v1

        # Calculate local distance difference 
        if kwargs["force_norm"]:
            # Normalize LDD by number of neighbors
            return np.linalg.norm(dv) / len(i1)
        else:
            return np.linalg.norm(dv)


    def calculate_ldd(self):
        """Calculate LDD"""
        self.ldd = self._calculate_deformation(self._calculate_ldd_residue)


    def _calculate_nd_residue(self, neigh_tensor1, neigh_tensor2, **kwargs):
        """Calculate neighborhood distance given a pair of neighborhood tensors"""
        # Rotate neighbourhood tensor and calculate Euclidean distance
        nd = np.linalg.norm(rotate_points(neigh_tensor2, neigh_tensor1) - neigh_tensor1)
        if kwargs["force_norm"]:
            return nd / len(i1)
        else:
            return nd


    def calculate_neighborhood_dist(self):
        """Calculate neighborhood_distance"""
        self.neighbor_distance = self._calculate_deformation(self._calculate_nd_residue)


    def _calculate_shear_residue(self, neigh_tensor1, neigh_tensor2, **kwargs):
        """Calculate shear strain given a pair of neighborhood tensors"""
        u1 = neigh_tensor1
        u2 = neigh_tensor2
        try:
            duu = u2 @ u2.T - u1 @ u1.T
            uu = np.linalg.inv(u1.T @ u1) 
            C = 0.5 * (uu @ u1.T @ duu @ u1 @ uu)
            return 0.5 * np.sum(np.diag(C@C) - np.diag(C)**2)
        except Exception as e:
            print(e)
            return np.nan


    def calculate_shear(self):
        """Calculate shear strain"""
        self.shear = self._calculate_deformation(self._calculate_shear_residue)


    def _calculate_strain_residue(self, neigh_tensor1, neigh_tensor2, **kwargs):
        """Calculate effective strain given a pair of neighborhood tensors"""
        # Rotate neighbourhood tensor and calculate Euclidean distance
        nt3 = rotate_points(neigh_tensor2, neigh_tensor1)
        if not kwargs["force_absolute"]:
            # Divide by length to get strain
            es = np.sum(np.linalg.norm(nt3 - neigh_tensor1, axis=1) / np.linalg.norm(neigh_tensor1, axis=1))
        else:
            es = np.sum(np.linalg.norm(nt3 - neigh_tensor1, axis=1))

        if not kwargs["force_nonorm"]:
            # Normalize ES by number of neighbors
            es /= len(neigh_tensor1)

        return es


    def calculate_strain(self):
        """Calculate effective strain"""
        self.strain = self._calculate_deformation(self._calculate_strain_residue)


    def _calculate_non_affine_residue(self, neigh_tensor1, neigh_tensor2, **kwargs):
        """Calculate non-affine strain given a pair of neighborhood tensors"""
        # Find the deformation gradient tensor, F
        F, residuals = np.linalg.lstsq(neigh_tensor1, neigh_tensor2, rcond=None)[:2]
        if not kwargs["force_nonorm"]:
            return np.nansum(residuals)
        else:
            return np.nansum(residuals) / len(neigh_tensor1)


    def calculate_non_affine(self):
        """Calculate non-affine strain"""
        self.non_affine = self._calculate_deformation(self._calculate_non_affine_residue)


    #######################################################
    ### Root-Mean-Square Deviation
    ### Superimpose structures and calculate RMSD
    def calculate_rmsd(self):
        """
        Calculate RMSD

        Cannot calculate RMSD with averaged neighborhoods,
        so if an AverageProtein object is passsed, we simply
        take the first Protein object and use it to calculate RMSD.
        """
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

        c1 = c1 - np.mean(c1, axis=0).reshape(1, 3)
        c2 = c2 - np.mean(c2, axis=0).reshape(1, 3)
        c3 = rotate_points(c2, c1)

        sd = np.sum((c1 - c3)**2, axis=1)

        self.rmsd_per_residue = np.zeros(len(idx)) * np.nan
        self.rmsd_per_residue[idx] = np.sqrt(sd)

        self.rmsd = np.sqrt(np.mean(sd))


