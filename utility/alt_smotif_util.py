import utility.io_util as io


def get_lowest_noe_energy(tasks):
    """
    :param tasks:
    :return:
    """
    noe_energy = []
    for entry in tasks:
        noe_energy.append(entry[5][3])
    half_len = (len(noe_energy) / 2.0)
    noe_energy = noe_energy[0:int(half_len)]
    return sum(noe_energy) / half_len


def compute_jobs(tasks):

    # print tasks[0][1]
    # print tasks[0][8]
    # print tasks[0][9]
    # Alt Smotif ['Alt_smotif', [1, 2, 'right']] number 9 in the list

    alt_smotif = tasks[0][9][1]

    if alt_smotif[-1] == 'right':
        alt_sse = alt_smotif[-2]
    else:
        alt_sse = alt_smotif[0]

    ss_profile = io.getSSprofilesFile()
    alt_sse_profile = ss_profile[alt_sse]
    # print alt_sse, alt_sse_profile
    jobs = []
    for i in range(0, len(tasks)):
        for j in range(0, len(alt_sse_profile)):
            jobs.append([i, j])
    # print jobs
    return jobs, alt_sse_profile


def getfromDB(pair, sse_ordered, database_cutoff):
    from utility.smotif_util import getSmotif, readSmotifDatabase
    s1 = sse_ordered[pair[0]]
    s2 = sse_ordered[pair[1]]
    s1_len = s1[1]
    s2_len = s2[1]
    smotif = getSmotif(s1, s2)
    return readSmotifDatabase(smotif, database_cutoff)


def getSmotifDB(sse_ordered, ss_profile, alt_smotif_log, pair, cutoff):

    if alt_smotif_log[-1] == 'right':
        sse_ordered[-1] = ss_profile
    else:
        sse_ordered[0] = ss_profile
    print pair, sse_ordered
    return getfromDB(pair, sse_ordered, cutoff)


def delete_last_sse(sse_coors, alt_smotif_log):

    if alt_smotif_log[-1] == 'right':
        return sse_coors[:-1]
    else:
        return sse_coors[1:]

