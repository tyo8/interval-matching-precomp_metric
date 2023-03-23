import ast
import argparse
import numpy as np


def get_prevalence_scores(verbose_match_list):
    prevalence_scores = [ np.mean(verbose_match["affinity"]) for verbose_match in verbose_match_list]
    return prevalence_scores

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Given the (collated) verbose matches from a set of bootstraps, compute and write per-cycle prevalence scores"
    )
    parser.add_argument(
        "-i", "--inpath", type=str, help="input filepath to (string literal of) verbose match list"
    )
    parser.add_argument(
        "-o", "--outpath", type=str, default=None, help="output filepath to prevalence scores"
    )
    parser.add_argument(
        "-v", "--verbose", default=False, action='store_true', help="Send a string literal of prevalence scores to stdout"
    )
    args = parser.parse_args()
    verbose_match_list_fname = args.inpath

    with open(verbose_match_list_fname,'r') as fin:
        verbose_match_list = ast.literal_eval(fin.read())

    prevalence_scores = get_prevalence_scores(verbose_match_list)

    if args.outpath:
        with open(args.outpath,'w') as fout:
            [fout.write(str(score)+'\n') for score in prevalence_scores]
    else:
        args.verbose=True

    if args.verbose:
        print("Prevalence scores dervied from " + args.inpath)
        print(prevalence_scores)
        print('')
