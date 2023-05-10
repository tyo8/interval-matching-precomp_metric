import gc
import os
import compute
import extract
import datetime
import numpy as np
from find_match import find_match
from benchmark_matching import benchmark_out
from create_matrices_image import create_matrices_image


def bootstrap_matching(
        dX=[], X=[], subidx_list=[], dim = 2, affinity_method = 'D', do_hom0 = True,
        verbose_out=False, verbose_figs = False, benchmark = True
        ):

    print('Verbose output set to: ' + str(verbose_out))
    print('')

    if benchmark:
        start = datetime.datetime.now()

    nb_X = dX.shape[0]
    out_X = compute.bars_tightreps(inp=dX, filename='dX', maxdim=dim, verbose = verbose_out)
    bars_X, reps_X, tight_reps_X, indices_X = extract.xtr_bars_reps_indices(out_X, do_hom0 = do_hom0, verbose = verbose_out)

        #_saveout(out_X, fname='phom_X.txt')    # osbsolete until relevant troubleshooting; phom_X should ofc be invaraint across bootstrappings, so we can recoup the write time

    list_matched_X_Y = {}
    list_affinity_X_Y = {}

    list_bars_reps_Y = []

    for i, subidx in enumerate(subidx_list) :
        # gives likely-unique identifier ('i'_'microsec') to bootstrapped intermediate saveouts to prevent cross-worker mixups
        stamp = str(i) + '_' + str(datetime.datetime.now())[-6:]
        print('Matching X to X_{} ...'.format(stamp))


        dZ,dY = _subsamp_dZ(dX, subidx)

        # Compute the lower distance matrices
        ldm_file_X, ldm_file_Y, ldm_file_Z, threshold = create_matrices_image(
                pairwise_Z=dZ, nb_X=nb_X, return_thr=True, write_mode=True,
                filename_X = 'dX', filename_Y = 'dX_'+stamp, filename_Z = 'dZ_'+stamp
                )

        # Image persistence - apply thresholding 
        out_X_Z = compute.image_bars(filename_X = ldm_file_X, filename_Z = ldm_file_Z, threshold = threshold, verbose = verbose_out)
        if benchmark:
            msg = 'Ripser-image (X-->Z) in iteration ' + str(i) + ' completed.'
            benchmark_out(msg, datetime.datetime.now() - start)
            
        bars_X_Z, indices_X_Z = extract.xtr_bars_indices(out_X_Z, do_hom0 = do_hom0, verbose = verbose_out)
        if benchmark:
            msg = 'Extraction of Ripser-image (X-->Z) in iteration ' + str(i) + ' completed.'
            benchmark_out(msg, datetime.datetime.now() - start)
            

        _saveout(out_X_Z, fname='phom_XZ_' + stamp + '.txt')

        out_Y_Z = compute.image_bars(filename_X = ldm_file_Y, filename_Z = ldm_file_Z, threshold = threshold, verbose = verbose_out)
        bars_Y_Z, indices_Y_Z = extract.xtr_bars_indices(out_Y_Z, do_hom0 = do_hom0, verbose = verbose_out)

        _saveout(out_Y_Z, fname='phom_YZ_' + stamp + '.txt')

        out_Y = compute.bars_tightreps(inp=dY, filename='dY_'+stamp, maxdim=dim, verbose = verbose_out)
        if benchmark:
            msg = 'Last Ripser (Y) in iteration ' + str(i) + ' completed.'
            benchmark_out(msg, datetime.datetime.now() - start)
            
        bars_Y, reps_Y, tight_reps_Y, indices_Y = extract.xtr_bars_reps_indices(out_Y, do_hom0 = do_hom0, verbose = verbose_out)
        if benchmark:
            msg = 'Extraction of last Ripser (Y) in iteration ' + str(i) + ' completed.'
            benchmark_out(msg, datetime.datetime.now() - start)
            

        _saveout(out_Y, fname='phom_Y_' + stamp + '.txt')

        list_bars_reps_Y += [ [bars_Y, reps_Y, tight_reps_Y] ]

        matched_X_Y, affinity_X_Y = find_match(bars_X, bars_X_Z, indices_X, indices_X_Z,
                                       bars_Y, bars_Y_Z, indices_Y, indices_Y_Z, dim = 1,
                                       affinity_method = affinity_method, check_Morse = False,
                                       check_ambiguous_deaths = False)
        if benchmark:
            msg = 'Cycle matching in iteration ' + str(i) + ' completed. Iteration over'
            benchmark_out(msg, datetime.datetime.now() - start)
            print('\n\n')

        list_matched_X_Y[i] = matched_X_Y
        list_affinity_X_Y[i] = affinity_X_Y

        if verbose_figs:
            print('sorry, no verbose figs yet!')
            ## have to figure out sibling/parent imports to make this work
            # from match.plot import show_matches
            # show_matches(X, list_Y[i], matched_X_Y[i], affinity_X_Y[i], tight_reps_X, tight_reps_Y, dim = dim, show_together = True)

    list_matched_X_Y = list(list_matched_X_Y.values())
    list_affinity_X_Y = list(list_affinity_X_Y.values())
    bars_reps_X = [bars_X, reps_X, tight_reps_X]
    
    return list_matched_X_Y, list_affinity_X_Y, bars_reps_X, list_bars_reps_Y

def _subsamp_dZ(dX, subsamp_idx):
    dXY = dX[:, subsamp_idx]
    dY = dXY[subsamp_idx, :]

    dZ = np.block([[dX, dXY], [dXY.T, dY]])
    return dZ, dY

def _saveout(output, fname='out.txt'):
    print('Saving Ripser output to: ' + os.path.abspath(fname))
    print('')
    with open(fname,'w') as fout:
        fout.write('\n'.join(output))


def multiple_matching(X, list_Y, dim = 1, verbose_figs = False, affinity_method = 'A') :

    out_X = compute.bars_tightreps(X)
    bars_X, reps_X, tight_reps_X, indices_X = extract.xtr_bars_reps_indices(out_X, do_hom0 = True)

    list_matched_X_Y = {}
    list_affinity_X_Y = {}

    list_bars_reps_Y = []

    for y, Y in enumerate(list_Y) :
        print('Matching X to Y_{} ...'.format(y))

        # Compute the lower distance matrices
        ldm_file_X, ldm_file_Y, ldm_file_Z, threshold = \
            create_matrices_image(X, Y, filename_X = 'X', filename_Y = 'Y', return_thr = True)

        # Image persistence - apply thresholding 
        out_X_Z = compute.image_bars(filename_X = ldm_file_X, filename_Z = ldm_file_Z, threshold = threshold)
        bars_X_Z, indices_X_Z = extract.xtr_bars_indices(out_X_Z, do_hom0 = True)

        out_Y_Z = compute.image_bars(filename_X = ldm_file_Y, filename_Z = ldm_file_Z, threshold = threshold)
        bars_Y_Z, indices_Y_Z = extract.xtr_bars_indices(out_Y_Z, do_hom0 = True)

        out_Y = compute.bars_tightreps(Y)
        bars_Y, reps_Y, tight_reps_Y, indices_Y = extract.xtr_bars_reps_indices(out_Y, do_hom0 = True)

        list_bars_reps_Y += [ [bars_Y, reps_Y, tight_reps_Y] ]

        matched_X_Y, affinity_X_Y = find_match(bars_X, bars_X_Z, indices_X, indices_X_Z,
                                       bars_Y, bars_Y_Z, indices_Y, indices_Y_Z, dim = 1,
                                       affinity_method = affinity_method, check_Morse = False,
                                       check_ambiguous_deaths = False)

        list_matched_X_Y[y] = matched_X_Y
        list_affinity_X_Y[y] = affinity_X_Y

        if verbose_figs :
            print('sorry, no verbose figs yet!')
            ## have to figure out sibling/parent imports to make this work
            # from match.plot import show_matches
            # show_matches(X, list_Y[y], matched_X_Y[y], affinity_X_Y[y], tight_reps_X, tight_reps_Y, dim = dim, show_together = True)

    list_matched_X_Y = list(list_matched_X_Y.values())
    list_affinity_X_Y = list(list_affinity_X_Y.values())
    bars_reps_X = [bars_X, reps_X, tight_reps_X]

    return list_matched_X_Y, list_affinity_X_Y, bars_reps_X, list_bars_reps_Y


def cross_matching(list_X, dim = 1, verbose_figs = False, affinity_method = 'A') :

    for i in range(len(list_X)) :
        for j in range(i) :
            aa = list_matched_X_Y[j,i]
            if len(aa) > 0 :
                list_matched_X_Y[i,j] = np.array(aa)[:,::-1].tolist() # reverse column order
            else :
                list_matched_X_Y[i,j] = []
            list_affinity_X_Y[i,j] = list_affinity_X_Y[j,i]
    return list_matched_X_Y, list_affinity_X_Y, list_bars_reps_X

##### SOME FUNCTIONS TO ENABLE TRACKING CYCLES


##### SOME FUNCTIONS TO ENABLE TRACKING CYCLES

def track_cycles_from_slice(list_matched_X_Y, list_affinity_X_Y, cycle, list_indices, initial_slice = 0):
    ''' From a list of matched cycles between consecutive slices, obtains a list in which we store the matches
        that track a particular cycle from a some chosen slice
        output = [[cycle, a],[a, b], [b, c] ...]
        rmk: set of indices does not include the initial slice, counts from the second slice studied'''

    tracked_cycle = []
    tracked_affinity = []
    current_match = []

    # initialise
    for i, match in enumerate(list_matched_X_Y[initial_slice]):
        if match[0] == cycle:
            tracked_cycle += [match]
            current_match = match
            tracked_affinity.append(list_affinity_X_Y[initial_slice][i])


    # track the cycle
    for i in list_indices:
        next_cycle = current_match[1]
        tracked_copy = tracked_cycle.copy()
        for j, match in enumerate(list_matched_X_Y[i]):
            if match[0] == next_cycle:
                #print(match)
                tracked_cycle += [match]
                current_match = match
                tracked_affinity.append(list_affinity_X_Y[i][j])
                continue
        if len(tracked_copy) == len(tracked_cycle): # in case we don't find a next cycle we stop tracking

            break

    return tracked_cycle, tracked_affinity
