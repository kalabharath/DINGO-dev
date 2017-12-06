import os

import io_util as io
import utility.smotif_util as sm

def enum(*sequential, **named):
    """

    :param sequential:
    :param named:
    :return:
    """
    # fake an enumerated type in Python

    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)


def getPairSSProfiles(s1, s2, ss_profile):
    """
    return corresponding +/- 2 profile
    definitons for the given given index
    order
    """
    s1_l = ss_profile[s1]
    s2_l = ss_profile[s2]
    return s1_l, s2_l


def getRunSeq():
    """
    generate run seq, a seq list of pairs of
    indexes of profiles for job scheduling
    """
    # TODO needed to generalize based on input data type
    map_route = []
    if os.path.isfile("ss_profiles.pickle"):
        ss_profiles = io.readPickle("ss_profiles.pickle")
    else:
        return False
    if os.path.isfile("contacts_route.pickle"):
        map_route = io.readPickle("contacts_route.pickle")
    elif os.path.isfile("pcs_route.pickle"):
        map_route = io.readPickle("pcs_route.pickle")
    elif os.path.isfile("rdc_route.pickle"):
        map_route = io.readPickle("rdc_route.pickle")
    print map_route
    s1, s2 = map_route[0][0], map_route[0][1]
    s1_list, s2_list = getPairSSProfiles(s1, s2, ss_profiles)

    run_seq = []
    for i in range(len(s1_list)):
        for j in range(len(s2_list)):
            run_seq.append([i, j])
    return run_seq


def getSSlist():
    map_route = []
    ss_profiles = io.readPickle("ss_profiles.pickle")
    if os.path.isfile("contacts_route.pickle"):
        map_route = io.readPickle("contacts_route.pickle")
    elif os.path.isfile("pcs_route.pickle"):
        map_route = io.readPickle("pcs_route.pickle")
    elif os.path.isfile("rdc_route.pickle"):
        map_route = io.readPickle("rdc_route.pickle")

    if map_route:
        s1, s2 = map_route[0][0], map_route[0][1]
        s1_list, s2_list = getPairSSProfiles(s1, s2, ss_profiles)
        return s1_list, s2_list
    else:
        return False, False


def getSSdef(index_array):
    """

    :param index_array:
    :return:
    """
    s1_list, s2_list = getSSlist()
    return s1_list[index_array[0]], s2_list[index_array[1]]


def getfromDB(previous_smotif, current_ss, direction, database_cutoff, stage):
    """
    from stage2
    :param database_cutoff:
    :param previous_smotif:
    :param current_ss:
    :param direction:
    :return:
    """

    if stage == 2:
        psmotif = []
        for entry in previous_smotif:
            if 'smotif_def' == entry[0]:
                psmotif = entry[-1]
        if direction == 'left':
            previous_ss = psmotif[0]
        else:
            previous_ss = psmotif[1]
    else:
        searched_smotifs = []
        for entry in previous_smotif:
            if 'smotif_def' == entry[0]:
                searched_smotifs = entry[-1]

        # ['smotif_def', [['helix', 6, 7, 5, 145, 150], ['helix', 23, 5, 1, 156, 178], ['strand', 5, 7, 8, 133, 137]]]

        if direction == 'left':
            previous_ss = searched_smotifs[0]
        else:
            previous_ss = searched_smotifs[-1]

    # current_ss, previous_ss
    if direction == 'left':  # double check this implementation
        smotif_def = sm.getSmotif(current_ss, previous_ss)
    else:
        smotif_def = sm.getSmotif(previous_ss, current_ss)

    return sm.readSmotifDatabase(smotif_def, database_cutoff), smotif_def


def generate_refinement_order(sse_array):

    import itertools
    from utility.smotif_util import array2string
    computed_pairs = []
    indices = list(itertools.combinations(range(len(sse_array)), 2))
    refine_pairs = []
    for pair in indices:
        if abs(pair[0] - pair[1]) == 1:
            pass
        else:
            t_array = array2string([sse_array[pair[0]], sse_array[pair[1]]])
            computed_pairs.append(t_array)
            refine_pairs.append(pair)

    return refine_pairs, computed_pairs

def generate_refinement_order2(sse_array, computed_pairs):

    import itertools
    from utility.smotif_util import array2string

    indices = list(itertools.combinations(range(len(sse_array)), 2))
    refine_pairs = []
    for pair in indices:
        if abs(pair[0] - pair[1]) == 1:
            pass
        else:
            t_array = array2string([sse_array[pair[0]], sse_array[pair[1]]])

            if t_array in computed_pairs:
                pass
            else:
                computed_pairs.append(t_array)
                refine_pairs.append(pair)

    return refine_pairs, computed_pairs





def orderSSE(previous_smotif, current_sse, direction, stage):
    """
    :param previous_smotif:
    :param current_sse:
    :param direction:
    :return:
    """
    if stage == 2:
        previous_seq = []
        for entry in previous_smotif:
            if 'smotif_def' == entry[0]:
                previous_seq = entry

        ordered_SSE = []
        for sse_array in previous_seq[1]:
            ordered_SSE.append(sse_array)

        if direction == 'left':
            ordered_SSE.insert(0, current_sse)
        else:
            ordered_SSE.append(current_sse)

        refine_pairs, computed_pairs = generate_refinement_order(ordered_SSE)
        return ordered_SSE, refine_pairs, computed_pairs
    else:
        computed_pairs = []
        for entry in previous_smotif:
            if 'Refine_Smotifs' == entry[0]:
                computed_pairs = entry[2]

        for entry in previous_smotif:
            # ['qcp_rmsd', transformed_coos, sse_ordered, rmsd]
            if 'qcp_rmsd' == entry[0]:
                previous_sse = entry[2]
                if direction == 'left':
                    previous_sse.insert(0, current_sse)
                else:
                    previous_sse.append(current_sse)

                refine_pairs, computed_pairs = generate_refinement_order2(previous_sse, computed_pairs)

                return previous_sse, refine_pairs, computed_pairs


def fetchNOEdata(previous_smotif):
    return previous_smotif[5][-3], previous_smotif[5][-2], previous_smotif[5][-1]