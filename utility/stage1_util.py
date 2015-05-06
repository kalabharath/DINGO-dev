import io_util as io


def enum(*sequential, **named):
    """

    :param sequential:
    :param named:
    :return:
    """
    #fake an enumerated type in Python

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
    ss_profiles = io.readPickle("ss_profiles.pickle")
    map_route = io.readPickle("contact_route.pickle")
    s1, s2 = map_route[0][0], map_route[0][1]
    s1_list, s2_list = getPairSSProfiles(s1, s2, ss_profiles)

    run_seq = []
    for i in range(len(s1_list)):
        for j in range(len(s2_list)):
            run_seq.append([i, j])
    return run_seq


def getSSlist():
    ss_profiles = io.readPickle("ss_profiles.pickle")
    map_route = io.readPickle("contact_route.pickle")
    s1, s2 = map_route[0][0], map_route[0][1]
    s1_list, s2_list = getPairSSProfiles(s1, s2, ss_profiles)
    return s1_list, s2_list

    # getRunSeq()