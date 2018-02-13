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
    smotif_def = getSmotif(s1, s2)
    return readSmotifDatabase(smotif_def, database_cutoff), sse_ordered, smotif_def


def getSmotifDB(sse_ordered, ss_profile, alt_smotif_log, pair, cutoff):

    if alt_smotif_log[-1] == 'right':
        sse_ordered[-1] = ss_profile
    else:
        sse_ordered[0] = ss_profile
    return getfromDB(pair, sse_ordered, cutoff)


def delete_last_sse(sse_coors, alt_smotif_log):

    if alt_smotif_log[-1] == 'right':
        return sse_coors[:-1]
    else:
        return sse_coors[1:]


def getSeq(coor_array, sse_ordered, aa_seq):

    one_letter = {'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'GLU': 'E', 'GLN': 'Q',
                  'ASP': 'D', 'ASN': 'N', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F', 'TYR': 'Y',
                  'ARG': 'R', 'LYS': 'K', 'SER': 'S', 'THR': 'T', 'MET': 'M', 'ALA': 'A',
                  'GLY': 'G', 'PRO': 'P', 'CYS': 'C', 'ASX': 'D', 'GLX': 'G', 'UNK': 'A'}
    concat_seq = ''
    for frag in coor_array:
        atom_num = 1
        for i in range(atom_num, len(frag[0]), 5):
            res = (frag[5][i])
            concat_seq = concat_seq+one_letter[res]

    native_sse_seq = ''
    for sse in sse_ordered:
        sse_seq = aa_seq[sse[4] - 1: sse[5]]
        native_sse_seq = native_sse_seq + sse_seq
    k = 0.0
    if len(concat_seq) != len(native_sse_seq):
        print "Something is wrong with extracting sequence information"
    for i in range(0, len(concat_seq)):

        if native_sse_seq[i] == concat_seq[i]:
            k += 1
    seq_id = (k / float(len(concat_seq))) * 100
    return concat_seq, seq_id


def altRMSDcutoff(smotif_def):
    from utility.smotif_util import getRMSDcutoff
    return getRMSDcutoff(smotif_def)
