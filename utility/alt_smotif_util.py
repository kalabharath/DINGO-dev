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
    elif alt_smotif[-1] == 'left':
        alt_sse = alt_smotif[0]
    else:
        print "Something is wrong with this logic"
        exit()

    ss_profile = io.readPickle("./ss_profiles.pickle")
    alt_sse_profile = ss_profile[alt_sse]
    # print alt_sse, alt_sse_profile
    jobs = []
    for i in range(0, len(tasks)):
        for j in range(0, len(alt_sse_profile)):
            jobs.append([i, j])
    # print jobs
    return jobs


