import re
import os
import sys
import warnings
import argparse
import numpy as np

sys.path.append("/scratch/tyoeasley/brain_representations/src_py")
from generate_subindex import tag_to_subidx


def summarize_diagram(phom_fname, img_flag=True, do_hom0=True, debug=False, write=False):
    with open(phom_fname,'r') as fin:
        phom_str = fin.read().split('\n')

    if img_flag:
        bars, indices = xtr_bars_indices(phom_str, do_hom0 = do_hom0, debug = debug)
        if write:
            _write_out(bars,'bars',phom_fname)
            _write_out(indices,'indices',phom_fname)
        else:
            return bars, indices

    else:
        bars, reps, tight_reps, indices = xtr_bars_reps_indices(phom_str, do_hom0 = do_hom0, debug = debug)
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



def xtr_bars_reps(out, do_hom0 = False, debug = True, maxdim=2) :
    '''This function converts the output of compute_bars_tight_reps into list of bars and reps, organised by dimension'''

    PH_nameline = _get_PH_nameline(out, maxdim=maxdim)

    bars = _init_dimdict(maxdim)
    reps = _init_dimdict(maxdim)
    tight_reps = _init_dimdict(maxdim)

    # reps 0 [ [[0],[1]], [[1],[2]], etc ]
    # reps 1 [ [ [0,1], [1,2], [2,3], [3,0] ], same]

    if do_hom0 :

        # 0-dim PH bars and reps
        dim = 0
        i = PH_nameline[dim]+1
        while i < PH_nameline[dim + 1] :
            bar, inf_flag = _get_bar(out[i],i)
            bars[dim] += [ bar ]
            rep = _get_h0rep(out[i], inf_flag)
            reps[dim] += [ rep ]
            i+=1

    # >=1-dim PH bars and reps
    for dim in range(1,maxdim+1):
        i = PH_nameline[dim]+1

        if debug:
            ### debug code ###
            print('')
            print('In: xtr_bars_reps')
            print(f'current line-read number={i}, dim={dim}')
            print(f"PH_nameline dictionary: \n{PH_nameline}")
            ### debug code ###

        while i < PH_nameline[dim + 1] :
            if i == len(out) - 1 : # trivial string ''
                break
            bar,_ = _get_bar(out[i],i)
            # all "finite" bars (no missing second endpoint)
            bars[dim] += [ bar ] 
            i += 1 # next line for tight reps
            tight_reps[dim] += [ _get_genreps(out[i]) ]
            i += 1 # again, next line for reps
            reps[dim] += [ _get_genreps(out[i]) ]
            i += 1

    if debug :
        print('phom_bars:')
        print(bars)
        print('')
        print('phom_reps:')
        print(reps)
        print('')
        print('phom_tightreps:')
        print(tight_reps)
        print('')

    return bars, reps, tight_reps



def xtr_bars_reps_indices(out, do_hom0 = False, debug = False, maxdim=2):
    '''This function converts the output of compute_bars_tight_reps into list of bars, representatives and indices of the persistence pairs,
    organised by dimension. REMARK: you need to use the modified version or ripser_tight_representative_cycles'''
    PH_nameline = _get_PH_nameline(out, maxdim=maxdim)

    bars = _init_dimdict(maxdim)
    reps = _init_dimdict(maxdim)
    tight_reps = _init_dimdict(maxdim)
    indices = _init_dimdict(maxdim)

    # reps 0 [ [[0],[1]], [[1],[2]], etc ]
    # reps 1 [ [ [0,1], [1,2], [2,3], [3,0] ], same]

    if do_hom0 :

        # 0-dim PH bars and reps
        dim = 0
        i = PH_nameline[dim]+1

        if debug:
            ### debug code ###
            print('')
            print('In: xtr_bars_reps_indices')
            print(f'current line-read number={i}, dim={dim}, \noutline={out[i]}')
            print(f"PH_nameline dictionary: \n{PH_nameline}")
            # print('Output for homology dimension ' + str(dim) + ':')
            # print(*out[i:PH_nameline[dim+1]], sep='\n')
            ### debug code ###

        while i < PH_nameline[dim + 1] :
            bar, inf_flag = _get_bar(out[i],i)
            bars[dim] += [ bar ]
            rep = _get_h0rep(out[i], inf_flag, lineno=i)
            reps[dim] += [ rep ]
            i+=1

    # 1-dim PH bars and reps
    for dim in range(1,maxdim):

        i = PH_nameline[dim]+1

        if debug:
            ### debug code ###
            print('')
            print(f'In: xtr_bars_reps_indices')
            print(f'current line-read number={i}, dim={dim}')
            print(f"PH_nameline dictionary: \n{PH_nameline}")
            ### debug code ###

        while i < PH_nameline[dim + 1] :
            if i == len(out) - 1 : # trivial string ''
                break
            bar,_ = _get_bar(out[i],i)
            # all "finite" bars (no missing second endpoint)
            bars[dim] += [ bar ] 

            # indices
            indices[dim] += [ _get_idx(out[i]) ]

            i += 1 # next line for tight reps
            tight_reps[dim] += [ _get_genreps(out[i]) ]
            i += 1 # again, next line for reps
            reps[dim] += [ _get_genreps(out[i]) ]
            i += 1

    if debug :
        print(bars)
        print(reps)
        print(tight_reps)

    return bars, reps, tight_reps, indices


def xtr_bars(out, do_hom0 = False, debug = False, maxdim=2) :
    ''' This function converts the output of compute_image_bars into list of bars organised by dimension
    (simpler version than xtr_bars_reps, no reps for image-persistence)'''
    PH_nameline = _get_PH_nameline(out, maxdim=maxdim)

    bars = _init_dimdict(maxdim)

    if do_hom0 :

        # 0-dim PH bars
        dim = 0
        i = PH_nameline[dim]+1
        while i < PH_nameline[dim + 1] :
            bar,_ = _get_bar(out[i],i)
            bars[dim] += [ bar ]
            i+=1

    # 1-dim PH bars
    dim = 1
    i = PH_nameline[dim]+1

    if debug:
        ### debug code ###
        print('')
        print(f'In: xtr_bars')
        print(f'current line-read number={i}, dim={dim}')
        print(f"PH_nameline dictionary: \n{PH_nameline}")
        ### debug code ###

    while i < PH_nameline[dim + 1] :
        if i == len(out) - 1 : # trivial string ''
            break
        bar,_ = _get_bar(out[i],i)
        bars[dim] += [ bar ]
        i+=1

    if debug :
        print(bars)

    return bars


def xtr_bars_indices(out, do_hom0 = False, debug = False, maxdim=2) :
    ''' This function converts the output of compute_image_bars into list of bars and indices of the persistence pairs organised by dimension
    (simpler version than xtr_bars_reps_indices, no reps for image-persistence). REMARK: need to use the modified version of ripser-image!'''
    PH_nameline = _get_PH_nameline(out, maxdim=maxdim)

    bars = _init_dimdict(maxdim)
    indices = _init_dimdict(maxdim)

    if do_hom0 :

        # 0-dim PH bars
        dim = 0
        i = PH_nameline[dim]+1
        while i < PH_nameline[dim + 1] :
            bar,_ = _get_bar(out[i],i)
            bars[dim] += [ bar ]
            i+=1

    # 1-dim PH bars
    dim = 1
    i = PH_nameline[dim]+1

    if debug:
        ### debug code ###
        print('')
        print(f'In: xtr_bars_indices)')
        print(f'current line-read number={i}, dim={dim}')
        print(f"PH_nameline dictionary: \n{PH_nameline}")
        ### debug code ###

    while i < PH_nameline[dim + 1] :
        #bars
        if i == len(out) - 1 : # trivial string ''
            break
        bar,_ = _get_bar(out[i],i)
        bars[dim] += [ bar ]
        #indices
        indices[dim] += [ _get_idx(out[i]) ]

        i += 1

    if debug :
        print(bars)
        print(indices)

    return bars, indices


def _get_PH_nameline(out, maxdim=2, debug=True):
    # find after which line it starts enumerating intervals in dim 0,1,2
    PH_nameline = {dim: len(out) for dim in range(maxdim+2)}
    dim_prologue1 = 'persistent homology intervals in dim '
    dim_prologue2 = 'persistence intervals in dim '
    for i,line in enumerate(out) :
        if line.startswith(dim_prologue1) or line.startswith(dim_prologue2):
            dim = int(line.split(':')[0][-1])
            PH_nameline[dim] = i
            if debug:
                ### debugging code ###
                print(f"Found naming line for dimension {dim} at line {i}: \"{line}\"")
                print(f"Updated PH_nameline dictionary: {PH_nameline}")
                ### debugging code ###

    return PH_nameline

def _init_dimdict(maxdim):
    dimdict = {dim: [] for dim in range(maxdim+1)}
    return dimdict


#  pulls birth-death interval from single line of text output
def _get_bar(line_out, line_num):
    inf_flag=False

    interval_pattern = r"\[(\d*.\d*),(\d*.\d*)\)"
    x = re.search(interval_pattern, line_out)
    if not x:
        sci_note_pattern=r"\b-?[1-9](?:\.\d+)?[Ee][-+]?\d+\b"
        x = re.findall(sci_note_pattern, line_out)     # check for scientific notation
        if not x :
            warnings.warn(f"no intervals found in output line number {line_num}: \'{line_out}\'")
            bar = []
        else:
            if len(x) > 1:
                bar = [float(numstr) for numstr in x]
            else:
                # this isn't quite right -- should also parse line_out for cases where we have
                #   (1) sci-notation birth but not death, or 
                #   (2) nonzero birth with sci-notation death (tho this would be *very* near diagonal so might be ok to overlook)
                bar = [0.] + [float(numstr) for numstr in x]

    else:
        if x.group(2) != ' ': # finite bars first
            bar = [float(x.group(1)), float(x.group(2))]
        else:
            inf_flag = True
            bar = [float(x.group(1)), np.inf]

    return bar, inf_flag


# pulls representative from single line of text output (assuming output comes from 0-dim homology)
def _get_h0rep(line_out, inf_flag, lineno=None, debug=True):
    if inf_flag:
        rep_pattern = r"\{\[(\d+)\].*\}"
    else:
        rep_pattern = r"\{\[(\d+)\], \[(\d+)\]\}"

    y = re.search(rep_pattern, line_out)
    if debug:
        ### debugging code ###
        if y is None:
            print("Found no matches")
            print(f"match pattern: {rep_pattern}")
            print(f"output line number: {lineno}")
            print(f"output line: {line_out}\n")
        ### debugging code ###

    if inf_flag:
        h0rep = [ [int(y.group(1))] ]
    else:
        h0rep = [ [int(y.group(1)), int(y.group(2))] ]
    return h0rep


# pulls representatives from single line of text output (assuming output does not come from 0-dim homology)
def _get_genreps(line_out):
    y = re.findall(r"\[(\d+),(\d+)\] \(\d*.\d*\)", line_out)
    y = [list(elem) for elem in y]
    genreps  = [[int(e[0]), int(e[1])] for e in y]
    return genreps


# pulls indices from single line of text output
def _get_idx(line_out):
    z = re.search(r"indices: (\d*)-(\d*)", line_out)
    if not z : 
        z = re.search(r"indices: \d*-\d*", line_out)
        if not z :
            raise ValueError("no indices found --- are you using the modified version of ripser-tight-representative-cycles?")
    idx = [int(z.group(1)), int(z.group(2))]
    return idx

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Extract persistence bars, filtration(?) indices, and cycle representatives (when possible) from rawtext Ripser output"
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
        "-v", "--debug", default=False, action='store_true', help="True if debug output"
    )

    args = parser.parse_args()

    summarize_diagram(
            args.phom_fname,
            img_flag = args.img_flag,
            do_hom0 = args.do_hom0,
            debug = args.debug,
            write = args.write
            )
