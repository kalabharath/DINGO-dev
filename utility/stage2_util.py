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

def makeTopPickle(previous_smotif_index):

    hits = []
    regex=str(previous_smotif_index)+"_*.pickle"
    file_list = glob.glob(regex)
    for f in file_list:
        thits = io.readPickle(f)
        for thit in thits:
            hits.append(thit)
    """
     dump_log = [['smotif',['seq_filter', 'smotif_seq', 'seq_identity', "blosum62_score"],
                ['contacts_filter','no_of_contacts', '%_of_contacts_observed'],
                ['PCS_filter', 'tensor_fits']]]
    """
    print len(hits)

    for hit in hits:
        smotif = hit[0]
        seq_filter = hit[1]
        print seq_filter


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


