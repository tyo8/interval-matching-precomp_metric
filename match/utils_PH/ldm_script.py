import sys
import argparse
import numpy as np
import create_matrices_image as cmi
sys.path.append('/scratch/tyoeasley/brain_representations/src_py')
from generate_subindex import tag_to_subidx

##########################################################################################################################################
### test code ###
# the function below produces fixed results that will behave mildly inconsistent across tags;
# this is ok, because this function only exists to allow code to run smoothly on test data and is not present in any important analyses
def reduce_subidx(subidx, N_new):
    approx_prop = float(len(subidx))/max([float(i) for i in subidx])
    subidx_long = [idx for idx in subidx if idx < N_new]

    new_length = int(N_new*approx_prop)
    subidx_new = subidx_long[:new_length]
    return subidx_new
### test code ###
##########################################################################################################################################

# Computes and saves subsampled distance & image distance matrices given an initial distance matrix and an ASCII-encoded bit array
if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Make and save image & subfiltration matrices and subsample matrices given distance matrix and subsample tag"
    )
    parser.add_argument(
        "-x", "--dX_fname", type=str, help="path to distance matrix"
    )
    parser.add_argument(
        "-t", "--tag", type=str, help="subsampling bit array (as ASCII)"
    )
    parser.add_argument(
        "-z", "--dZ_fname", type=str, help="outpath to image distance matrix"
    )
    parser.add_argument(
        "-y", "--dY_fname", type=str, help="outpath to subsampled distance matrix"
    )
    parser.add_argument(
        "-i", "--dXZ_fname", type=str, help="outpath to embedded image distance matrix"
    )
    parser.add_argument(
        "-j", "--dYZ_fname", type=str, help="outpath to embedded subsampled image distance matrix"
    )
    args = parser.parse_args()

    dX = np.genfromtxt(args.dX_fname)
    nb_X = dX.shape[0]

    subidx = tag_to_subidx(args.tag)

    if len(subidx) >= nb_X:
        subidx = reduce_subidx(subidx,nb_X)

    dZ, dY = cmi._subsamp_dZ(dX, subidx)

    cmi._write_out(args.dY_fname, dY)
    cmi.create_matrices_image(
            pairwise_Z=dZ, 
            nb_X=nb_X,
            filename_X=args.dXZ_fname,
            filename_Y=args.dYZ_fname,
            filename_Z=args.dZ_fname
            )
