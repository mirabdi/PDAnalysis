import numpy as np

### "P" is the set of points to be mapped to "Q"
def rotate_points(P, Q):
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    V = Vt.T
    D = np.linalg.det(V @ U.T)
    E = np.array([[1, 0, 0], [0, 1, 0], [0, 0, D]])
    R = V @ E @ U.T
    Pnew = np.array([R @ p for p in P])
    return Pnew


### Returns a np.ndarray of positions where two sequences differ
def get_mutation_position(s1, s2):
    return np.where(np.array(list(s1)) != np.array(list(s2)))[0]


def get_shared_indices(idx1, idx2):
    i1 = [j for j, k in enumerate(idx1) if k in idx2]
    i2 = [j for j, k in enumerate(idx2) if k in idx1]
    return i1, i2


