import os
import re
import sys
import argparse
import numpy as np
import extract
sys.path.append('/scratch/tyoeasley/brain_representations/src_py')
from generate_subindex import tag_to_subidx


def summarize_diagram(phom_fname, img_flag=True, do_hom0=True, verbose=False, write=False):
    with open(phom_fname,'r') as fin:
        phom_str = fin.read().split('\n')

    if img_flag:
        bars, indices = extract.xtr_bars_indices(phom_str, do_hom0 = do_hom0, verbose = verbose)
        if write:
            _write_out(bars,'bars',phom_fname)
            _write_out(indices,'indices',phom_fname)
        else:
            return bars, indices

    else:
        bars, reps, tight_reps, indices = extract.xtr_bars_reps_indices(phom_str, do_hom0 = do_hom0, verbose = verbose)
        if write:
            _write_out(bars,'bars',phom_fname)
            _write_out(reps,'reps',phom_fname)
            _write_out(tight_reps,'tight_reps',phom_fname)
            _write_out(indices,'indices',phom_fname)
        else:
            return bars, reps, tight_reps, indices


def _write_out(var, varlabel, phom_fpath):
    savedir   = os.path.dirname(phom_fpath)
    phom_fname= os.path.basename(phom_fpath)
    var_fname = phom_fname.replace('phom', varlabel)
    var_fpath = os.path.join(savedir, var_fname)
    # print(var_fpath)
    with open(var_fpath,'w') as fout:
        varstring = var.__str__().replace('inf','2e308')     # replace 'inf' with a true literal (may break in future versions of Python)
        fout.write(varstring)


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Make and save image & subfiltration matrices and subsample matrices given distance matrix and subsample tag"
    )
    parser.add_argument(
        "-x", "--phom_fname", type=str, help="filepath to distance matrix"
    )
    parser.add_argument(
        "-0", "--do_hom0", default=False, action='store_true', help="True if output is expected to contain H0 data"
    )
    parser.add_argument(
        "-w", "--write", default=False, action='store_true', help="True if derived variables are written to file; if False, derived variables are returned"
    )
    parser.add_argument(
        "-i", "--img_flag", default=False, action='store_true', help="True if parsing image persistence output"
    )
    parser.add_argument(
        "-v", "--verbose", default=False, action='store_true', help="True if verbose output"
    )

    args = parser.parse_args()

    summarize_diagram(
            args.phom_fname,
            img_flag = args.img_flag,
            do_hom0 = args.do_hom0,
            verbose = args.verbose,
            write = args.write
            )
