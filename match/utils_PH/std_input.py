import numpy as np

def phom(X=[], Y=[], pairwise_Z=[], nb_X=[]):
    if not np.any(pairwise_Z):
        assert np.any(X), "If pairwise_Z is not provided, point clouds X and Y must both be given"
        assert np.any(Y), "If pairwise_Z is not provided, point clouds X and Y must both be given"
        Z = np.vstack((X,Y))
        nb_X = len(X)
        pairwise_Z = np.sqrt( np.sum( (Z[:, None, :] - Z[None, :, :])**2, axis=-1) )
    else:
        assert nb_X, "If pairwise_Z is given and X and Y are not, nb_X must be specified."

    return X, Y, pairwise_Z, nb_X

