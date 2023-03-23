import compute
import numpy as np

def duplicates_list(values) : # duplicates_list([1,2,3,1,2,2]) = [1, 2, 2]
    seen = set()
    dupes = [val for val in values if val in seen or seen.add(val)]
    return dupes

def find_occurences_list(values, val) :
    indices = [i for i, x in enumerate(values) if x == val]
    return indices


def find_match(bars_X, bars_X_Z, indices_X, indices_X_Z, bars_Y, bars_Y_Z, indices_Y, indices_Y_Z, dim = 1, affinity_method = 'A',
               check_Morse = False, check_ambiguous_deaths = False) :
    ''' This funtion find the matches between the barcodes of X and Y providing the barcodes of their image-persistence modules in the union.
    Affinity score is automatically set to A but can be changed. Optiona outputs to check if the filtrations provided are Morse and if 
    there are image-bars sharing death times in the barcodes.'''

    # initialize the match and affinity variables as dictionaries indexed by the indices of bars_XZ
    matched_X_Y = {idx: [] for idx, bar in enumerate(bars_X[dim])}
    affinity_X_Y = {idx: [] for idx, bar in enumerate(bars_X[dim])}
    verbose_matches = [{
        "dim":dim,
        "barX": bar,
        "barY": [],
        "deathZ": [],
        "matchedXY":[],
        "affinity":[0]} for bar in bars_X[dim]]

    # consider all image-bars
    births_X_Z = [bar[0] for bar in bars_X_Z[dim]]
    births_Y_Z = [bar[0] for bar in bars_Y_Z[dim]]
    deaths_X_Z = [bar[1] for bar in bars_X_Z[dim]]
    deaths_Y_Z = [bar[1] for bar in bars_Y_Z[dim]]
    ## unzip

    # consider normal bars
    births_X = [bar[0] for bar in bars_X[dim]]
    births_Y = [bar[0] for bar in bars_Y[dim]]
    deaths_X = [bar[1] for bar in bars_X[dim]]
    deaths_Y = [bar[1] for bar in bars_Y[dim]]
    ## unzip

    if check_Morse :
        # adding noise to your point clouds does not solve the following exceptions.
        # It will create a distance matrix with unique values, so that only adding 1 edge at a time 
        # but possibly many triangles that kill cycles simultaneously in Rips complexes
        if len(duplicates_list(deaths_X_Z)) > 0 :
            print('Found duplicate deaths in X_Z')
        if len(duplicates_list(deaths_Y_Z)) > 0 :
            print('Found duplicate deaths in Y_Z')
        if len(duplicates_list(births_X)) > 0 :
            print('Found duplicate births in X') # should never happen for unique distance values
        if len(duplicates_list(births_Y)) > 0 :
            print('Found duplicate births in Y') # should never happen for unique distance values


    considered_deaths_X_Z = set(deaths_X_Z)
    considered_deaths_Y_Z = set(deaths_Y_Z)

    # find common (considered) deaths in image
    common_deaths = considered_deaths_X_Z.intersection(considered_deaths_Y_Z)


    if check_ambiguous_deaths :
        if set(duplicates_list(deaths_X_Z)).intersection(set(duplicates_list(deaths_Y_Z))) != set() :
            print('Found common duplicate deaths in X_Z and Y_Z!!!')
        if set(duplicates_list(deaths_X_Z)).intersection(common_deaths) != set() :
            print('Found duplicate death in X_Z common with Y_Z')
        if set(duplicates_list(deaths_Y_Z)).intersection(common_deaths) != set() :
            print('Found duplicate death in Y_Z common with X_Z')

    # determine ambiguous deaths
    ambiguous_deaths_X_Z = set(duplicates_list(deaths_X_Z)).intersection(common_deaths)
    ambiguous_deaths_Y_Z = set(duplicates_list(deaths_Y_Z)).intersection(common_deaths)
    ambiguous_deaths = ambiguous_deaths_X_Z.union(ambiguous_deaths_Y_Z)
    if ambiguous_deaths != set() :
        print('We will solve ambiguous deaths matching.')

    # now, find common births with the normal bars

    # First case: non-ambiguous matching
    # in this case, deaths in X_Z and Y_Z are unique so matching can be made without ambiguity (even if duplicate deaths in X or in Y)

    for death in common_deaths.difference(ambiguous_deaths) :
        oXZ = deaths_X_Z.index(death)
        oYZ = deaths_Y_Z.index(death)
        birth_X = births_X_Z[oXZ]
        birth_Y = births_Y_Z[oYZ]

        # Now we match with the persistence bars of X and Y
        Occ_X = find_occurences_list(births_X, birth_X)
        Occ_Y = find_occurences_list(births_Y, birth_Y)

        if len(Occ_X) == 1 and len(Occ_Y) == 1: # if there are no ambiguous births
            a = births_X.index(birth_X) # unique 
            b = births_Y.index(birth_Y) # unique 

            matched_X_Y[a].append(b)

            if b:
                affinity_val = compute.affinity(birth_X, deaths_X[a], death, birth_Y, deaths_Y[b], affinity_method = affinity_method)
            else:
                affinity_val = -1
                b = 0

            affinity_X_Y[a].append(affinity_val)
            verbose_matches = _update_verbose(verbose_matches, birth_X, deaths_X, death, birth_Y, deaths_Y, a, b, affinity_val)
        else:
            # the way we are computing persistent homology, the indices of the 
            # persistent homology of Y can be compared with the indices in the
            # image - persistent homology of Y inside Z

            for oX in Occ_X:
                pos_index_X = indices_X[dim][oX][0]
                pos_index_XZ = indices_X_Z[dim][oXZ][0]
                if pos_index_X == pos_index_XZ:
                    a = oX
            if 'a' not in locals():
                print("oXZ = " + str(oXZ))
                print("Occ_X = " + str(Occ_Xz))
                print("pos_index_X = " + str(pos_index_X))
                print("pos_index_XZ = " + str(pos_index_XZ))
                raise Exception("birth index 'a' in X homology is unset -- this should be impossible.")

            for oY in Occ_Y:
                pos_index_Y = indices_Y[dim][oY][0]
                pos_index_YZ = indices_Y_Z[dim][oYZ][0] # the bars are presented in the same order, independently on how we arrange X an Y
                if pos_index_Y == pos_index_YZ:
                    b = oY
            if 'b' not in locals():
                print(f"\nWARNING: match index \'b\' unset in dimension {dim} at oYZ={oYZ}")
                b = ''

            matched_X_Y[a].append(b)

            if b:
                affinity_val = compute.affinity(birth_X, deaths_X[a], death, birth_Y, deaths_Y[b], affinity_method = affinity_method)
            else:
                affinity_val = -1
                b = 0

            affinity_X_Y[a].append(affinity_val)
            verbose_matches = _update_verbose(verbose_matches, birth_X, deaths_X, death, birth_Y, deaths_Y, a, b, affinity_val)

    # Second case: ambiguous matching

    for death in ambiguous_deaths :
        # detect the indices of the bars with ambiguous death times
        Occ_XZ = find_occurences_list(deaths_X_Z, death)
        Occ_YZ = find_occurences_list(deaths_Y_Z, death)

        for i, oXZ in enumerate(Occ_XZ) :
            # extract negative index of the image persistence bar of X
            neg_index_XZ = indices_X_Z[dim][oXZ][1]
            for j, oYZ in enumerate(Occ_YZ) :
                # extract negative index of the image persistence bar of Y
                neg_index_YZ = indices_Y_Z[dim][oYZ][1]
                if neg_index_XZ == neg_index_YZ:
                    # match the image bars when the indices coincide
                    birth_X = births_X_Z[oXZ] # ! not i
                    birth_Y = births_Y_Z[oYZ] # ! not j

                    # Now we match with the persistence bars of X and Y
                    Occ_X = find_occurences_list(births_X, birth_X)
                    Occ_Y = find_occurences_list(births_Y, birth_Y)

                    if len(Occ_X) == 1 and len(Occ_Y) == 1: # if there are no ambiguous births
                        a = births_X.index(birth_X) # unique 
                        b = births_Y.index(birth_Y) # unique 
                        matched_X_Y[a].append(b)

                        if b:
                            affinity_val = compute.affinity(birth_X, deaths_X[a], death, birth_Y, deaths_Y[b], affinity_method = affinity_method)
                        else:
                            affinity_val = -1
                            b = 0

                        affinity_X_Y[a].append(affinity_val)
                        verbose_matches = _update_verbose(verbose_matches, birth_X, deaths_X, death, birth_Y, deaths_Y, a, b, affinity_val)
                    else:
                        # the way we are computing persistent homology, the indices of the 
                        # persistent homology of Y can be compared with the indices in the
                        # image - persistent homology of Y inside Z
                        for oX in Occ_X:
                            pos_index_X = indices_X[dim][oX][0]
                            pos_index_XZ = indices_X_Z[dim][oXZ][0]
                            if pos_index_X == pos_index_XZ:
                                a = oX
                        if 'a' not in locals():
                            print("oXZ = " + str(oXZ))
                            print("Occ_X = " + str(Occ_Xz))
                            print("pos_index_X = " + str(pos_index_X))
                            print("pos_index_XZ = " + str(pos_index_XZ))
                            raise Exception("birth index 'a' in X homology is unset -- this should be impossible.")

                        for oY in Occ_Y:
                            pos_index_Y = indices_Y[dim][oY][0]
                            pos_index_YZ = indices_Y_Z[dim][oYZ][0] # the bars are presented in the same order, independently on how we arrange X an Y
                            if pos_index_Y == pos_index_YZ:
                                b = oY
                        if 'b' not in locals():
                            print(f"\nWARNING: match index \'b\' unset in dimension {dim} at oYZ={oYZ}")
                            b = ''

                        matched_X_Y[a].append(b)

                        if b:
                            affinity_val = compute.affinity(birth_X, deaths_X[a], death, birth_Y, deaths_Y[b], affinity_method = affinity_method)
                        else:
                            affinity_val = -1
                            b = 0 

                        affinity_X_Y[a].append(affinity_val)
                        verbose_matches = _update_verbose(verbose_matches, birth_X, deaths_X, death, birth_Y, deaths_Y, a, b, affinity_val)


    return matched_X_Y, affinity_X_Y, verbose_matches

def _update_verbose(verbose_match_list, birth_X, deaths_X, deathZ, birth_Y, deaths_Y, a, b, aff_val):
    entry = verbose_match_list[a]

    barX_new = [birth_X, deaths_X[a]]
    assert np.allclose(np.asarray(entry["barX"]), np.asarray(barX_new)), "Error: conflict between barX entries " + str(entry["barX"]) + " and " + str(barX_new)

    entry["barY"].append([birth_Y, deaths_Y[b]])
    entry["deathZ"].append(deathZ)
    entry["matchedXY"].append([a,b])

    if entry["affinity"] == [0]:
        entry["affinity"] = []
    entry["affinity"].append(aff_val)

    return verbose_match_list
