import glob, os

import io_util as io


def enum(*sequential, **named):
    """

    :param sequential:
    :param named:
    :return:
    """
    # fake an enumerated type in Python

    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)


def checkFile(i):
    regex = str(i) + "_*.pickle"
    file_list = glob.glob(regex)
    if len(file_list) > 0:
        return True
    else:
        return False


def getNextSmotif(map_route):
    # [[6, 7, 'start'], [5, 6, 'left'], [4, 5, 'left'], [3, 4, 'left'], [2, 3, 'left'], [1, 2, 'left'], [0, 1, 'left']
    for i in range(0, len(map_route)):
        if not checkFile(i):
            return i, map_route[i]


def scoreCombination(score_list):
    """

    :param score_list:
    :return:
    """

    import itertools
    min_score = 999
    combi_list = list(itertools.combinations(score_list, 2))

    for combi in combi_list:
        c1 = combi[0]
        c2 = combi[1]
        if c1 and c2:
            if c1+c2 < min_score:
                min_score = c1 + c2
    return min_score

def scoreCombination4t(score_list):
    """

    :param score_list:
    :return:
    """

    import itertools
    min_score = 999
    combi_list = list(itertools.combinations(score_list, 3))

    for combi in combi_list:
        c1 = combi[0]
        c2 = combi[1]
        c3 = combi[2]
        if c1 and c2 and c3:
            if c1+c2+c3 < min_score:
                min_score = c1 + c2 + c3
    return min_score

def getNchiSum(pcs_filter):

    """

    :param pcs_filter:
    :return:
    """

    snchi = 999.999

    tensors = pcs_filter[1]
    #print len(tensors)

    if len(tensors) == 2:
        snchi = 0
        for tensor in tensors:
            # tag = tensor[0]
            nchi = tensor[1]
            # axrh = tensor[2]
            snchi += nchi

    if len(tensors) == 3:
        # Scoring three tags, get lowest Nchi for 2
        score_list = []
        for tensor in tensors:
            score_list.append(tensor[1])
        snchi = scoreCombination(score_list)

    if len(tensors) >= 4:
        # For 4 tags, get lowest Nchi for 3
        score_list = []
        for tensor in tensors:
            score_list.append(tensor[1])
        snchi = scoreCombination4t(score_list)

    if len(tensors) == 1:
        # Discourage single tag scoring by returning high score
        return snchi

    return snchi


def makeTopPickle(previous_smotif_index, num_hits):
    hits = []
    regex = str(previous_smotif_index) + "_*_*.pickle"
    file_list = glob.glob(regex)
    #print file_list
    for f in file_list:
        thits = io.readPickle(f)
        for thit in thits:
            hits.append(thit)
    """
    identifiers: smotif, smotif_def, seq_filter, contacts_filter, PCS_filter,

    """

    new_dict = {}

    for hit in hits:
        for entry in hit:
            if entry[0] == 'smotif':
                name = entry[1][0]
            if entry[0] == 'seq_filter':
                seq_filter = entry
                smotif_seq = seq_filter[1]
            if entry[0] == 'contacts_filter':
                contacts_filter = entry
            if entry[0] == 'PCS_filter':
                pcs_data = entry
                Nchi = getNchiSum(pcs_data)
                new_dict.setdefault(Nchi, []).append(hit)

    keys = new_dict.keys()
    keys.sort()

    non_redundant = {}
    seqs = []
    for i in range(0, len(keys)):
        for entry in new_dict[keys[i]][0]:
            if entry[0] == 'smotif':
                name = entry[1][0]
            if entry[0] == 'seq_filter':
                seq_filter = entry
                smotif_seq = seq_filter[1]
            if entry[0] == 'contacts_filter':
                contacts_filter = entry
            if entry[0] == 'PCS_filter':
                pcs_data = entry
                Nchi = getNchiSum(pcs_data)
        if smotif_seq not in seqs:
            seqs.append(smotif_seq)
            print name, Nchi
            non_redundant.setdefault(Nchi, []).append(new_dict[keys[i]][0])

    keys = non_redundant.keys()
    keys.sort()
    dump_pickle = []
    try:
        for i in range(0, num_hits):
            dump_pickle.append(non_redundant[keys[i]])
            print "final sele", non_redundant[keys[i]][0][0][1][0]
    except:
        print "Could only extract ", i
        num_hits = i
    io.dumpPickle(str(previous_smotif_index) + "_tophits.pickle", dump_pickle)
    #delete_old = "rm "+str(previous_smotif_index)+"_*_*.pickle"
    #os.system(delete_old)
    print "actual number in top hits ", len(dump_pickle)
    return range(num_hits)


def getRunSeq(num_hits):
    """
    generate run seq, a seq list of pairs of
    indexes of profiles for job scheduling
    """

    ss_profiles = io.readPickle("ss_profiles.pickle")
    map_route = io.readPickle("pcs_route.pickle")
    next_index, next_smotif = getNextSmotif(map_route)

    direction = next_smotif[-1]
    if direction == 'left':
        next_ss_list = ss_profiles[next_smotif[0]]
    else:
        next_ss_list = ss_profiles[next_smotif[1]]
    # get and make a list of top 10(n) of the previous run
    top_hits = makeTopPickle(next_index - 1, num_hits)  # send the previous Smotif index
    if top_hits:
        run_seq = []
        for i in range(len(top_hits)):
            for j in range(len(next_ss_list)):
                run_seq.append([i, j])
        return run_seq, next_index

def getPreviousSmotif(index):
    map_route = io.readPickle("pcs_route.pickle")
    next_index, next_smotif = getNextSmotif(map_route)
    top_hits = io.readPickle(str(next_index - 1) + "_tophits.pickle")  # Read in previous index hits
    # print len(top_hits)
    return top_hits[index][0]

def getSS2(index):
    ss_profiles = io.readPickle("ss_profiles.pickle")
    # map_route = io.readPickle("contact_route.pickle")
    map_route = io.readPickle("pcs_route.pickle")
    next_index, next_smotif = getNextSmotif(map_route)
    direction = next_smotif[-1]
    if direction == 'left':
        next_ss_list = ss_profiles[next_smotif[0]]
    else:
        next_ss_list = ss_profiles[next_smotif[1]]

    return next_ss_list[index], direction


def rename_pickle(index):
    import glob, os

    file_list = glob.glob("tx_*.pickle")
    for file in file_list:
        mv_cmd = "mv " + file + " " + str(index) + file[2:]
        os.system(mv_cmd)
    return True
