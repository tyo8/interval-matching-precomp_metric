import os
import re
import ast
import argparse
from find_match import find_match


def find_match_wrapper(
        read=False, phomX_fname='', phomY_fname='', phomXZ_fname='', phomYZ_fname='',
        barsX=None, barsY=None, barsXZ=None, barsYZ=None, 
        idxX=None, idxY=None, idxXZ=None, idxYZ=None,
        dim=1, affinity_method='C', check_Morse=False, check_ambiguous_deaths=False
        ):
    if read:
        assert phomX_fname, 'Must supply filepath (type) of phomX data when in read mode'
        assert phomY_fname, 'Must supply filepath (type) of phomY data when in read mode'
        assert phomXZ_fname,'Must supply filepath (type) of phomXZ data when in read mode'
        assert phomYZ_fname,'Must supply filepath (type) of phomYZ data when in read mode'

        barsX= _read_in(phomX_fname, 'bars')
        barsY= _read_in(phomY_fname, 'bars')
        barsXZ=_read_in(phomXZ_fname,'bars')
        barsYZ=_read_in(phomYZ_fname,'bars')

        idxX= _read_in(phomX_fname,  'indices')
        idxY= _read_in(phomY_fname,  'indices')
        idxXZ=_read_in(phomXZ_fname, 'indices')
        idxYZ=_read_in(phomYZ_fname, 'indices')
    else:
        assert barsX, 'Must supply barsX data when not in read mode'
        assert barsY, 'Must supply barsY data when not in read mode'
        assert barsXZ,'Must supply barsXZ data when not in read mode'
        assert barsYZ,'Must supply barsYZ data when not in read mode'
        assert idxX,  'Must supply idxX data when not in read mode'
        assert idxY,  'Must supply idxY data when not in read mode'
        assert idxXZ, 'Must supply idxXZ data when not in read mode'
        assert idxYZ, 'Must supply idxYZ data when not in read mode'

    matched, affinity, verbose_matches = find_match(
        barsX, barsXZ, idxX, idxXZ,
        barsY, barsYZ, idxY, idxYZ,
        dim=dim, affinity_method=affinity_method, 
        check_Morse=check_Morse, check_ambiguous_deaths=check_ambiguous_deaths
            )

    return matched, affinity, verbose_matches



def _pull_tag(phom_fpath):
    phom_fname = os.path.basename(phom_fpath).split('.')[0]
    tag = re.sub('phom[X-Z]?[X-Z]?_','',phom_fname)
    return tag

def _read_in(phom_fpath, varlabel):
    savedir   = os.path.dirname(phom_fpath)
    phom_fname= os.path.basename(phom_fpath)
    var_fname = phom_fname.replace('phom', varlabel)
    var_fpath = os.path.join(savedir, var_fname)
    with open(var_fpath,'r') as fin:
        var = ast.literal_eval(fin.read())
    return var

def _write_out(var, varlabel, phom_fpath):
    savedir   = os.path.dirname(phom_fpath)
    phom_fname= os.path.basename(phom_fpath)
    var_fname = phom_fname.replace('phom', varlabel)
    var_fpath = os.path.join(savedir, var_fname)
    print(var_fpath)
    with open(var_fpath,'w') as fout:
        varstring = var.__str__().replace('inf','2e308')     # replace 'inf' with a true literal (may break in future versions of Python)
        fout.write(varstring)

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Make and save image & subfiltration matrices and subsample matrices given distance matrix and subsample tag"
    )
    parser.add_argument(
        "-x", "--phomX_fname", type=str, help="filepath to X (fixed) persistent homology"
    )
    parser.add_argument(
        "-y", "--phomY_fname", type=str, help="filepath to Y (subsampled) persistent homology"
    )
    parser.add_argument(
        "-D", "--dim", default=1, type=int, help="dimension in which to compute affinity"
    )
    parser.add_argument(
        "-o", "--outdir", default="", type=str, help="output directory for matched indices and affinity scores"
    )
    parser.add_argument(
        "-a", "--affinity_method", default='C', type=str, help="selects from various cycle pairings from which to compute the Jaccard index-based cycle affinity measure"
    )
    parser.add_argument(
        "-M", "--Morse", default=False, action='store_true', help="flag to check Morse condition"
    )
    parser.add_argument(
        "-d", "--check_deaths", default=False, action='store_true', help="flag to check for ambiguous death values during cycle-matching"
    )
    args = parser.parse_args()

    if args.outdir:
        outdir = args.outdir
    else:
        outdir = os.path.join(os.path.dirname(args.phomX_fname), 'matching')

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    
    phomXZ_fname = args.phomY_fname.replace('phomY','phomXZ')
    phomYZ_fname = args.phomY_fname.replace('phomY','phomYZ')

    matched, affinity, verbose_matches = find_match_wrapper(read=True,
            phomX_fname=args.phomX_fname, phomY_fname=args.phomY_fname, 
            phomXZ_fname=phomXZ_fname, phomYZ_fname=phomYZ_fname,
            dim=args.dim, affinity_method=args.affinity_method, 
            check_Morse=args.Morse, check_ambiguous_deaths=args.check_deaths
            )

    print("Saving matching data to: ")
    tag = _pull_tag(args.phomY_fname)
    dummy_filename = os.path.join(outdir,'phom_dim' + str(args.dim) + '_' + tag + '.txt')
    _write_out(matched, 'matched', dummy_filename)
    _write_out(affinity,'affinity',dummy_filename)
    _write_out(verbose_matches, 'verbose_match', dummy_filename)
