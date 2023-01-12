import compute
import create_matrices_image
import find_match
import std_input
import extract


def matching(X=[], Y=[], pairwise_Z=[], nb_X=[], maxdim = 2, verbose_figs = False, affinity_method = 'A', check_Morse = False) :

    X, Y, pairwise_Z, nb_X = std_input.phom(X=X, Y=Y, pairwise_Z=pairwise_Z, nb_X=nb_X)
    # Compute the lower distance matrices
    # NOTE: change this so we hold the ldm variables in memory instead of read/writing
    ldm_fname_X, ldm_fname_Y, ldm_fname_Z, threshold = create_matrices_image(
            X=X, Y=Y, pairwise_Z=dZ, nb_X=nb_X, return_thr = True,
            write_mode=True, filename_X = 'X', filename_Y = 'Y'
            )

    # Image persistence - apply thresholding (bug fixed)
    out_X_Z = compute.compute_image_bars(
            filename_X = ldm_fname_X, filename_Z = ldm_fname_Z, threshold = threshold, maxdim=maxdim,
            )
    bars_X_Z, indices_X_Z = extract.extract_bars_indices(out_X_Z, maxdim=maxdim, do_hom0 = False)

    out_Y_Z = compute.compute_image_bars(
            filename_X = ldm_fname_Y, filename_Z = ldm_fname_Z, threshold = threshold, maxdim=maxdim,
            )
    bars_Y_Z, indices_Y_Z = extract.extract_bars_indices(out_Y_Z, maxdim=maxdim, do_hom0 = False)

    # Persistent homology
    out_X = compute.compute_bars_tightreps(inp = None, filename = ldm_fname_X, maxdim=maxdim)
    # this way we obtain the same bars but the vertices are indexes in the same way as in the image persistence
    bars_X, reps_X, tight_reps_X, indices_X = extract.extract_bars_reps_indices(out_X, do_hom0 = False, maxdim=maxdim)

    out_Y = compute.compute_bars_tightreps(inp = None, filename = ldm_fname_Y, maxdim=maxdim)
    bars_Y, reps_Y, tight_reps_Y, indices_Y = extract.extract_bars_reps_indices(out_Y, do_hom0 = False, maxdim=maxdim)

    matched_X_Y, affinity_X_Y = find_match(bars_X, bars_X_Z, indices_X, indices_X_Z,
                                           bars_Y, bars_Y_Z, indices_Y, indices_Y_Z, dim = 1,
                                           affinity_method = affinity_method, check_Morse = check_Morse,
                                           check_ambiguous_deaths = False)

    if verbose_figs :
        ## have to figure out sibling/parent imports to make this work
        # from match.plot import show_matches
        # show_matches(X,Y,matched_X_Y, affinity_X_Y, tight_reps_X, tight_reps_Y, dim = dim, show_together = True)

    return matched_X_Y, affinity_X_Y, (bars_X, reps_X, tight_reps_X), (bars_Y, reps_Y, tight_reps_Y)


