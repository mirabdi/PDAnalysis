from collections import defaultdict
from pathlib import Path

from Bio import Seq, SeqIO, SeqRecord, Align
from Bio.Data import IUPACData
from Bio.PDB import PDBParser, MMCIF2Dict
import numpy as np
import pandas as pd


#############################################################
### Load coordinates from PDB file


### Read 3-letter code in title case
def parse_3letter(x):
    """Convert 3-letter amino acid code to a single-letter"""
    return IUPACData.protein_letters_3to1.get(x[0].upper() + x[1:].lower(), 'X')


def parse_pdb_coordinates(path, chain='', model=0, all_atom=False):
    """Parse coordinates, residue indices, sequence, and bfactor/plddt from a PDB file"""
    parser = PDBParser(QUIET=True)
    chains = list(parser.get_structure('', path)[model])
    coord, idx, seq, bfac = [], [], [], []

    for ch in chains:
        if (ch.get_id() == chain) or (chain == ''):
            for residue in ch:
                try:
                    h, i, _ = residue.get_id()
                    # Ignore HETATM entries
                    if h.strip() != '':
                        continue
                    aa = residue.resname
                    if all_atom:
                        for atom in residue:
                            coord.append(atom.coord)
                            idx.append(i)
                            seq.append(parse_3letter(aa))
                            bfac.append(atom.bfactor)
                    else:
                        if 'CA' not in residue:
                            continue
                        ca = residue['CA'].coord
                        coord.append(ca)
                        idx.append(i)
                        seq.append(parse_3letter(aa))
                        bfac.append(residue['CA'].bfactor)
                except Exception as e:
                    print(f"{path}\n{e}")
                    continue
    return [np.array(x) for x in [coord, idx, seq, bfac]]




#############################################################
### Load coordinates from mmCIF file


def parse_mmcif_coordinates(path, chain=''):
    """Parse coordinates, residue indices, sequence, and bfactor/plddt from a mmCIF file"""
    mmcif = reformat_mmcif_dict(MMCIF2Dict.MMCIF2Dict(path))
    df = pd.DataFrame(data=mmcif['_atom_site'])

    # Extract CA atoms
    # model '1'
    # ignore entires with incorrect sequence IDs
    # if there are multiple configurations, only choose 'A'
    df = df.loc[(df.label_atom_id=='CA')&(df.pdbx_PDB_model_num=='1')&(df.label_seq_id!='.')&(df.label_alt_id!='B')]

    # Extract specific CHAIN if arg is provided
    if len(chain):
        df = df.loc[(df.label_asym_id==chain)]

    coord = np.array([df[x].values for x in ['Cartn_x', 'Cartn_y', 'Cartn_z']]).T.astype(float)
    idx = df['label_seq_id'].values.astype(int)
    seq = np.array([parse_3letter(aa) for aa in df.label_comp_id])
    bfac = df["B_iso_or_equiv"].values.astype(float)
    return [coord, idx, seq, bfac]


def reformat_mmcif_dict(mmcif):
    """Reformat mmCIF dictionary to a nested dictionary"""
    new_dict = defaultdict(dict)
    for k, v in mmcif.items():
        if '.' not in k:
            new_dict[k] = v 
            continue
        k1 = k.split('.')[0]
        k2 = '.'.join(k.split('.')[1:])
        new_dict[k1][k2] = v 
    return new_dict



#############################################################
### Load SEQRES 

def load_pdb_seqres(path, chain=''):
    """Load SEQRES sequence from PDB file"""
    for record in SeqIO.parse(path, "pdb-seqres"):
        if not len(chain):
            return str(record.seq)
        elif record.annotations['chain'] == chain:
            return str(record.seq)


def load_mmcif_seqres(path, chain=''):
    """Load SEQRES sequence from mmCIF file"""
    mmcif = reformat_mmcif_dict(MMCIF2Dict.MMCIF2Dict(path))
    df = pd.DataFrame(data=mmcif['_pdbx_poly_seq_scheme'])

    # Extract specific CHAIN if arg is provided
    if len(chain):
        df = df.loc[(df.asym_id==chain)]

    seq = np.array([parse_3letter(aa) for aa in df.mon_id])
    return ''.join(seq)


#############################################################
### Load PDB file, reorder indices so they start from zero,
### and fill in missing coordinates / Bfactor with NaN values.


### Load SEQRES sequence, and use that to get a more complete
### description of the protein structure that preserves residues
### that are missing atoms.
### Doesn't work for some rare weird cases that you get in the PDB:
###     e.g. microheterogeneity??? (multiple amino acids for one site, somehow; eg. 1eis)
def load_and_fix_pdb_data(path, chain=''):
    """
    Read PDB / mmCIF file

    Load the SEQRES sequence, and match the atomic coordinates to
    indices based on the SEQRES sequence.
    """
    ext = Path(path).suffix

    # Load the sequence from the SEQRES part
    # Load the coords, sequence, etc., from the ATOM part
    if ext in ['.pdb', '.ent']:
        seqres = load_pdb_seqres(path, chain)
        xyz, idx, seq, bfac = parse_pdb_coordinates(path, chain)

    elif ext == '.cif':
        seqres = load_mmcif_seqres(path, chain)
        xyz, idx, seq, bfac = parse_mmcif_coordinates(path, chain)

    # If the ATOM sequence is longer than the SEQRES sequence,
    # then there is a problem
    if len(seqres) < len(seq):
        raise Exception("ERROR: the number of alpha-carbons exceeds the number of SEQRES residues!" + \
                        "\n\tPlease check your input files for errors.")

    # If there are no missing atoms, the two sequences will be equal.
    # In this case, the indices will be an integer series starting at 0
    if seqres == ''.join(seq):
        idx = np.arange(len(seq))
    else:
        # If there are missing atoms, try to align the SEQRES / ATOM sequences
        # to get the correct indices
        is_clear, idx = match_xyz_indices_to_seqres(seqres, xyz, seq)
        # If not "is_clear" (if any atom positions are ambiguous),
        # ignore ambiguous atoms
        if not is_clear:
            if len(idx):
                xyz, seq, bfac = [x[idx] for x in [xyz, seq, bfac]]
            else:
                raise Exception(f"Error reading pdb file\n\t{path}")

    # Fill in missing coordinates with nan values
    xyz_nan = np.zeros((len(seqres), 3), float) * np.nan
    xyz_nan[idx] = xyz 

    return [np.array(x) for x in [xyz, idx, list(seqres), bfac, xyz_nan]]


### Match SEQRES sequence to the sequence obtained from 
### the atomic coordinates. Do NOT allow any mismatches, only gaps.
def match_xyz_indices_to_seqres(seqres, xyz, seq):
    """Match the SEQRES sequence to the sequence obtained from ATOM entries"""
    # Find all residues that are connected along the backbone,
    # and cluster them into unbroken stretches of amino acids
    seq_clusters = find_neighbours(seq, xyz)
    # Run a strict alignment algorithm that discards candidate alignments
    # if they do not agree provide the same set of unbroken sequences
    # of amino acids (sequence clusters)
    candidates = align_sequences(seqres, ''.join(seq), seq_clusters)
    if len(candidates) > 1:
        # If there is more than one alignment, return a single index
        # that corresponds to the positions in the first alignment.
        return False, resolve_ambiguity(candidates)

    elif len(candidates) == 1:
        # If there is only one candidate...
        return True, np.where(np.array(list(candidates[0])) != '-')[0]

    else:
        # If there are no candidates identified with the strict alignment algorithm,
        # run without enforcing equivalence of sequence clusters
        candidates = align_sequences(seqres, ''.join(seq))
        if len(candidates) > 1:
            return False, resolve_ambiguity(candidates)

        elif len(candidates) == 1:
            return True, np.where(np.array(list(candidates[0])) != '-')[0]

        else:
            # If there are still no candidates, return False
            return False, []


### Distance between alpha-carbons in neighbouring
### amino acids ought to be about 3.8 AA.
### Cutoff is higher than this, due to inaccuracies in the PDB.
### See "12ca", positions 125-127 as an example
### Occasionally, this method fails because a string of disordered residues
### are missing, yet the amino acids bookending the missing residues are indexed beside
### each other; not actually that rare (e.g., 4s34_A)
def find_neighbours(seq, xyz, cut=4.3):
    """Find backbone neighbors based on their Alpha-carbon distances"""
    D = np.linalg.norm(xyz[:-1] - xyz[1:], axis=1)
    seq_clusters = []
    cluster = seq[0]
    for s, d in zip(seq[1:], D):
        if d < cut:
            cluster = cluster + s
        else:
            seq_clusters.append(cluster)
            cluster = s
    seq_clusters.append(cluster)
    return seq_clusters


def align_sequences(seqres, seq, seq_clusters=''):
    """Align SEQRES sequences to ATOM sequences"""
    candidates = []

    # Convert sequence clusters to set
    if not isinstance(seq_clusters, str):
        seqA = set(seq_clusters)

    # Loop through candidate alignments
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    # Set mismatch penalty to negative infinity,
    # since no mismatches are not allowed
    aligner.mismatch_score = -np.inf
    for align in aligner.align(seqres, seq):
        s1, s2 = align.sequences
        # If sequence clusters are provided, only allow candidates that
        # include all sequence clusters as unbroken sequences
        if not isinstance(seq_clusters, str):
            # Break alignment into sequence clusters by splitting
            # at gaps
            clusters = [c for c in s2.split('-') if len(c)]
            # If all of the previously identified sequence clusters
            # are found in the alignment, add the candidate to the output
            if seqA.issubset(set(clusters)):
                candidates.append(s2)
        else:
            candidates.append(s2)
    return candidates


### Discard any ambiguous residues.
# If there are ambiguous regions, they will be
# as long as the number of candidates, so we can
# account for the indices
# These cases are rare: the algorithm is only called for 3 PDB entries out of ~4000.
def resolve_ambiguity(cand):
    """Resolve ambiguous cases where there is more than one alignment between SEQRES and ATOM sequences"""
    cand = np.array([list(x) for x in cand])
    N = len(cand)
    idx = []
    icount = -1
    # Loop through sequence positions (including gaps)
    for i in range(cand.shape[1]):
        # If no gaps, there is no ambiguity; count the index
        if '-' not in cand[:,i]:
            icount += 1
            idx.append(icount)
        # If there are gaps, ignore the index if there is a gap
        # in the first candidate alignment
        else:
            if cand[0,i] != '-':
                icount += 1
#               idx.append(icount)
    return np.array(idx)




