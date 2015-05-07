import  io_util as io
import glob

def enum(*sequential, **named):
    """

    :param sequential:
    :param named:
    :return:
    """
    #fake an enumerated type in Python

    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

def checkFile(i):
    regex=str(i)+"_*.pickle"
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

def getNchiSum(pcs_filter):
    #['PCS_filter', [[0, 0.04950933584267763, [[29.943566759232581, 7.8711499604434305], [35.840777738719297, 8.0029596961619411]]], [1, 0.001111114208027107, [[142.57986684837724, 86.977485137698238], [91.67464082627454, 54.850648632962333]]]]]
    tensors = pcs_filter[1]
    snchi =0
    for tensor in tensors:
        tag = tensor[0]
        nchi = tensor[1]
        axrh = tensor[2]
        snchi += nchi
    return snchi

def makeTopPickle(previous_smotif_index):

    hits = []
    regex=str(previous_smotif_index)+"_*.pickle"
    file_list = glob.glob(regex)
    for f in file_list:
        thits = io.readPickle(f)
        for thit in thits:
            hits.append(thit)
    """
                    0                                  1
     dump_log = [['smotif',['seq_filter', 'smotif_seq', 'seq_identity', "blosum62_score"],
                                           2
                ['contacts_filter','no_of_contacts', '%_of_contacts_observed'],
                              3
                ['PCS_filter', 'tensor_fits']]]
    """
    print len(hits)
    new_dict={}
    seqs = []
    for hit in hits:
        smotif = hit[0]
        seq_filter = hit[1]
        smotif_seq = seq_filter[1]
        if smotif_seq not in seqs:
            seqs.append(smotif_seq)
            pcs_filter = hit[3]
            nchi_sum = getNchiSum(pcs_filter)
            new_dict.setdefault(nchi_sum, []).append(hit)
    keys = new_dict.keys()
    keys.sort()
    print len(keys)
    dump_pickle = []
    for i in range(0,10):
        dump_pickle.append(new_dict[keys[i]])
    io.dumpPickle(str(previous_smotif_index)+"_tophits.pickle", dump_pickle)
    return True

def getRunSeq():

    """
    generate run seq, a seq list of pairs of
    indexes of profiles for job scheduling
    """
    ss_profiles = io.readPickle("ss_profiles.pickle")
    map_route = io.readPickle("contact_route.pickle")
    next_index, next_smotif = getNextSmotif(map_route)

    direction = next_smotif[-1]
    if direction == 'left':
        next_ss_list = ss_profiles[next_smotif[0]]
    else:
        next_ss_list = ss_profiles[next_smotif[1]]


    ##get and make a list of top 10(n) of the previous run
    if makeTopPickle(next_index -1): # send the previous Smotif index

        return True


