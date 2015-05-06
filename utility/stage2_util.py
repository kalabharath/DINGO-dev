import  io_util as io


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
    import glob
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


def getRunSeq():

    """
    generate run seq, a seq list of pairs of
    indexes of profiles for job scheduling
    """
    ss_profiles = io.readPickle("ss_profiles.pickle")
    map_route = io.readPickle("contact_route.pickle")
    next_index, next_smotif = getNextSmotif(map_route)
    print next_index, next_smotif
    direction = next_smotif[-1]
    if direction == 'left':
        next_ss_list = ss_profiles[next_smotif[0]]
    else:
        next_ss_list = ss_profiles[next_smotif[1]]

    print next_ss_list
    ##get and make a list of top 10(n) of the previous run

    return next_ss_list
