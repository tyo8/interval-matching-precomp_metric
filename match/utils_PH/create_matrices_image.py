import os
import std_input
import numpy as np


def create_matrices_image(X=[], Y=[], pairwise_Z=[], nb_X=[], return_thr=False,
        write_mode=True, filename_X='X', filename_Y='Y', filename_Z='dZ'):
    '''Function to create the matrices for the computation of image-persistence so that we can compare the indices of the persistence pairs'''
    X, Y, pairwise_Z, nb_X = std_input.phom(X=X, Y=Y, pairwise_Z=pairwise_Z, nb_X=nb_X)

    maxi = np.max(pairwise_Z)
    pairwise_X = pairwise_Z.copy()
    pairwise_Y = pairwise_Z.copy()
    pairwise_X[nb_X:] = 2 * maxi + 1 # add min offset 1 because maxi can be small
    pairwise_Y[:,:nb_X] = 2 * maxi + 1 #observe here the change wrt the previous line
    threshold = 2 * maxi + 0.5 # to make sure we still include people <= maxi (in case maxi = 0... paranoia lol)

    # so that later we can apply thresholding in Ripser-image (bug fixed)
    # threshold = 2 * maxi # not 2 * maxi - 1 as maxi could be very small

    if write_mode:
        ldm_X = _write_out(filename_X, pairwise_X)
        ldm_Y = _write_out(filename_Y, pairwise_Y)
        ldm_Z = _write_out(filename_Z, pairwise_Z)

    else:
        ldm_X = pairwise_X
        ldm_Y = pairwise_Y
        ldm_Z = pairwise_Z

    if return_thr :
        return ldm_X, ldm_Y, ldm_Z, threshold

    return ldm_X, ldm_Y, ldm_Z


def _subsamp_dZ(dX, subsamp_idx):
    dXY = dX[:, subsamp_idx]
    dY = dXY[subsamp_idx, :]

    dZ = np.block([[dX, dXY], [dXY.T, dY]])
    return dZ, dY


def _write_out(fname, square_pdist):
    if not fname.endswith('.ldm'):
        fname = fname + '.ldm'

    ldm_fname = os.path.abspath(fname)
    with open(ldm_fname, "w") as f: # erase possibly pre-existing
        for i in range(len(square_pdist)):
            f.write(', '.join([ str(x) for x in square_pdist[i,:i]]))
            f.write('\n')

        f.close()

    return ldm_fname
