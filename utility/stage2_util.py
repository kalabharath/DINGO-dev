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
        #tag = tensor[0]
        nchi = tensor[1]
        #axrh = tensor[2]
        snchi += nchi
    return snchi

def makeTopPickle(previous_smotif_index, num_hits):

    hits = []
    regex=str(previous_smotif_index)+"_*_*.pickle"
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
                              3                        4
                ['PCS_filter', 'tensor_fits']['smotif_def',[ss1,ss2]]]
    """

    new_dict={}
    seqs = []
    for hit in hits:
        print hit
        seq_filter = hit[1]
        smotif_seq = seq_filter[1]
        contacts_filter = hit[2]
        if smotif_seq not in seqs:
            seqs.append(smotif_seq)
            dict_key_score = seq_filter[2]+contacts_filter[2]
            new_dict.setdefault(dict_key_score, []).append(hit)
    keys = new_dict.keys()
    keys.sort()
    dump_pickle = []
    for i in range(0,num_hits):
        dump_pickle.append(new_dict[keys[i]])
    io.dumpPickle(str(previous_smotif_index)+"_tophits.pickle", dump_pickle)
    return range(num_hits)

def getRunSeq(num_hits):

    print num_hits
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

    print next_ss_list
    ##get and make a list of top 10(n) of the previous run
    top_hits = makeTopPickle(next_index -1, num_hits) # send the previous Smotif index
    if top_hits:
        run_seq = []
        for i in range(len(top_hits)):
            for j in range(len(next_ss_list)):
                run_seq.append([i, j])
        return run_seq, next_index

def makeTopPickle3(previous_smotif_index, num_hits):

    hits = []
    regex=str(previous_smotif_index)+"_*_*.pickle"
    file_list = glob.glob(regex)
    for f in file_list:
        thits = io.readPickle(f)
        for thit in thits:
            hits.append(thit)
    """
                            0                                  1
    dump_log.append([transformed_coos,['seq_filter', csse_seq, seq_identity, blosum62_score],
                            2                                                           3
                ['contacts_filter', no_of_contacts, percent_of_satisfied_contacts], sse_ordered])

    """
    #print len(hits)
    new_dict={}
    seqs = []
    for hit in hits:
        #print hit
        seq_filter = hit[1]
        smotif_seq = seq_filter[1]
        contacts_filter = hit[2]
        if smotif_seq not in seqs:
            seqs.append(smotif_seq)
            dict_key_score = seq_filter[2]+contacts_filter[2]
            new_dict.setdefault(dict_key_score, []).append(hit)
    keys = new_dict.keys()
    keys.sort()
    dump_pickle = []
    try:
        for i in range(0,num_hits):
            dump_pickle.append(new_dict[keys[i]])
    except:
        print "could only make ", i
    io.dumpPickle(str(previous_smotif_index)+"_tophits.pickle", dump_pickle)
    return range(num_hits)



def getRunSeq3(num_hits):

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
    top_hits = makeTopPickle3(next_index -1, num_hits) # send the previous Smotif index
    if top_hits:
        run_seq = []
        for i in range(len(top_hits)):
            for j in range(len(next_ss_list)):
                run_seq.append([i, j])
        return run_seq, next_index

def getPreviousSmotif(index):
     map_route = io.readPickle("contact_route.pickle")
     next_index, next_smotif = getNextSmotif(map_route)
     top_hits = io.readPickle(str(next_index-1)+"_tophits.pickle") #Read in previous index hits
     #print len(top_hits)
     return top_hits[index][0]

def getPreviousSmotif3(index):
     map_route = io.readPickle("contact_route.pickle")
     next_index, next_smotif = getNextSmotif(map_route)
     top_hits = io.readPickle(str(next_index-1)+"_tophits.pickle") #Read in previous index hits
     #print len(top_hits)
     return top_hits[index][0], top_hits[index][-1]


def getSS2(index):
    ss_profiles = io.readPickle("ss_profiles.pickle")
    map_route = io.readPickle("contact_route.pickle")
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
        mv_cmd = "mv "+file+" "+str(index)+file[2:]
        os.system(mv_cmd)
    return True
