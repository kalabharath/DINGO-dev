import collections
import glob
import math
import os

import utility.io_util as io


def checkFile(i):
    # regex = str(i) + "_*.pickle"
    regex = str(i) + "_*.gzip"
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
            if c1 + c2 < min_score:
                min_score = c1 + c2
    return min_score


def scoreCombination4t(score_list):
    """

    :param score_list:
    :return:
    """

    import itertools
    min_score = 999.999
    combi_list = list(itertools.combinations(score_list, 3))

    for combi in combi_list:
        c1 = combi[0]
        c2 = combi[1]
        c3 = combi[2]
        if c1 and c2 and c3:
            if c1 + c2 + c3 < min_score:
                min_score = c1 + c2 + c3
    return min_score


def getNchiSum(pcs_filter, stage):
    """

    :param pcs_filter:
    :return:
    """
    snchi = 999.999
    tensors = pcs_filter[1]

    if len(tensors) == 1:
        # Discourage single tag scoring by returning high score
        return 999.999  # discourage double tag score only for 4 tags

    if len(tensors) == 2 and stage == 2:  # stage 2
        snchi = 0
        for tensor in tensors:
            nchi = tensor[1]
            snchi += nchi
            # return 999.999 #discourage double tag score only for 4 tags

    if len(tensors) == 3 and stage <= 3:  # stage 2 & 3
        # Scoring three tags, get lowest Nchi for 2
        score_list = []
        for tensor in tensors:
            score_list.append(tensor[1])
        snchi = scoreCombination(score_list)

    if len(tensors) >= 4 and stage <= 3:  # stage 2,3 & 4
        # For 4 tags, get lowest Nchi for 3
        score_list = []
        for tensor in tensors:
            score_list.append(tensor[1])
        snchi = scoreCombination4t(score_list)
        snchi /= 100.0  # artificially increase the priority

    if len(tensors) >= 4 and stage == 4:  # stage 4
        # For 4 tags, get lowest Nchi for 3
        score_list = []
        for tensor in tensors:
            score_list.append(tensor[1])
        snchi = scoreCombination4t(score_list)

    if len(tensors) == 3 and stage == 99:
        score_list = []
        for tensor in tensors:
            score_list.append(tensor[1])
        snchi = score_list[0] + score_list[1] + score_list[2]

    if stage == 999:
        if len(tensors) < 4:
            snchi = 999.999
        else:
            score_list = []
            for tensor in tensors:
                score_list.append(tensor[1])
            snchi = score_list[0] + score_list[1] + score_list[2] + score_list[3]

    return snchi


def rdcSumChi(rdc_data, stage):
    snchi = 999.999
    tensors = rdc_data[1]
    if len(tensors) == 1:
        for tensor in tensors:
            return tensor[0]

    if len(tensors) == 2 and stage == 2:
        snchi = 0
        for tensor in tensors:
            nchi = tensor[0]
            snchi += nchi
    if len(tensors) > 2:
        snchi = 0
        for tensor in tensors:
            nchi = tensor[0]
            snchi += nchi

    return snchi


def makeTopPickle(previous_smotif_index, num_hits, stage):
    """
    Concatenate data from all of the threads, organize, remove redundancies, rank
     and extract top hits as defined
    :param previous_smotif_index:
    :param num_hits:
    :param stage:
    :return:
    """
    hits = []
    # regex = str(previous_smotif_index) + "_*_*.pickle"
    regex = str(previous_smotif_index) + "_*_*.gzip"
    file_list = glob.glob(regex)
    for f in file_list:
        t_hits = io.readGzipPickle(f)
        for t_hit in t_hits:
            hits.append(t_hit)
    """
    identifiers: smotif, smotif_def, seq_filter, contacts_filter, PCS_filter, qcp_rmsd, Evofilter
                 RDC_filter, NOE_filter
    """

    new_dict = collections.defaultdict(list)
    rdc_constant = 0.0
    for hit in hits:
        # thread_data contains data from each search and filter thread.
        # initialize total score array
        total_score = {}
        for data_filter in hit:

            if data_filter[0] == 'PCS_filter':
                pcs_data = data_filter
                pcsscore = getNchiSum(pcs_data, stage)
                total_score['pcs_score'] = pcsscore

            if data_filter[0] == 'Ref_RMSD':
                total_score['rmsd_score'] = data_filter[1]

            if data_filter[0] == 'RDC_filter':
                rdc_data = data_filter
                #Nchi = rdcSumChi(rdc_data, stage)
                log_likelihood = data_filter[2]
                rdc_tensors = data_filter[1]
                for tensor in rdc_tensors:
                    rdc_constant = rdc_constant + tensor[0]
                rdc_constant = rdc_constant * 1e-10
                total_score['rdc_score'] = log_likelihood

            if data_filter[0] == 'NOE_filter':
                noe_probability = data_filter[1]
                log_likelihood = -1* (math.log(noe_probability))
                total_score['noe_score'] = log_likelihood

                # calculate the total score and append the hit
        if total_score:
            keys = total_score.keys()
            keys = ['noe_score','rdc_score']
            #keys = ['rmsd_score']
            tscore = 0
            for key in keys:
                tscore = tscore + total_score[key]
            tscore = tscore + rdc_constant
            if tscore < 999.999:
                new_dict[tscore].append(hit)

    # ************************************************
    # Exclude the redundant entries and rank top hits
    # ************************************************

    keys = new_dict.keys()
    keys.sort()

    # Exclude the redundant data.

    # non_redundant = {}
    non_redundant = collections.defaultdict(list)
    seqs = []
    smotif_seq = ''
    count_hits = 0
    for i in range(0, len(keys)):
        entries = new_dict[keys[i]]
        for entry in entries:
            for ent in entry:
                if ent[0] == 'seq_filter':
                    seq_filter = ent
                    smotif_seq = seq_filter[1]
            if smotif_seq not in seqs:
                seqs.append(smotif_seq)
                non_redundant[keys[i]].append(entry)
                count_hits += 1
        if count_hits >= num_hits:
            break

    # Rank top hits and dump the data
    keys = non_redundant.keys()
    keys.sort()

    dump_pickle = []
    print "Dumping data to disk"
    count_top_hits = 0
    while (True):
        for key in keys:
            entries = non_redundant[key]
            for entry in entries:
                dump_pickle.append(entry)
                # print "final selected Smotif: ", entry[0][1][0][0], "with score: ", key
                print "final sele", entry[0][1][0][0], key
                count_top_hits += 1
            if count_top_hits >= num_hits:
                break
        if count_top_hits >= num_hits:
            break
        else:
            print "could only extract ", count_top_hits
            break

    # io.dumpPickle(str(previous_smotif_index) + "_tophits.pickle", dump_pickle)
    io.dumpGzipPickle(str(previous_smotif_index) + "_tophits.gzip", dump_pickle)
    print "actual number in top hits ", len(dump_pickle)
    return range(count_top_hits)

def makeTopPickle2(previous_smotif_index, num_hits, stage):
    """
    Concatenate data from all of the threads, organize, remove redundancies, rank
     and extract top hits as defined
    :param previous_smotif_index:
    :param num_hits:
    :param stage:
    :return:
    """
    hits = []
    # regex = str(previous_smotif_index) + "_*_*.pickle"
    regex = str(previous_smotif_index) + "_*_*.gzip"
    file_list = glob.glob(regex)
    for f in file_list:
        t_hits = io.readGzipPickle(f)
        for t_hit in t_hits:
            hits.append(t_hit)
    """
    identifiers: smotif, smotif_def, seq_filter, contacts_filter, PCS_filter, qcp_rmsd, Evofilter
                 RDC_filter, NOE_filter
    """

    new_dict = collections.defaultdict(list)

    for hit in hits:
        # thread_data contains data from each search and filter thread.

        if hit[4][0] == 'NOE_filter':
            no_of_noes = hit[4][2]
            new_dict[no_of_noes].append(hit)

    keys = new_dict.keys()
    keys.sort()
    keys.reverse()
    # Rank based on NOE energy
    non_redundant = collections.defaultdict(list)
    reduced_dump_log = []
    seqs = []
    count_hits = 0
    for i in range(len(keys)):
        entries = new_dict[keys[i]]
        if len(entries) == 1:  # There is only one entry in this no_of_noes bin just check of existing sequences and move on
            smotif_seq = entries[0][3][1]
            if smotif_seq not in seqs:
                seqs.append(smotif_seq)
                reduced_dump_log.append(entries[0])
                print "final sele", entries[0][0][1][0][0], keys[i]
                count_hits += 1
        else:
            t_log = collections.defaultdict(list)
            for hit in entries:  # filter on noe_energy
                if hit[4][0] == 'NOE_filter':
                    noe_energy = hit[4][3]
                    noe_energy = round(noe_energy, 2)
                    t_log[noe_energy].append(hit)
            noe_energy_bins = t_log.keys()
            noe_energy_bins.sort()

            for j in range(len(noe_energy_bins)):  # filter on RDC score
                t2_log = collections.defaultdict(list)
                hits = t_log[noe_energy_bins[j]]
                for hit in hits:
                    if hit[5][0] == 'RDC_filter':
                        rdc_tensors = hit[5][1]
                        rdc_score = 0
                        for tensor in rdc_tensors:
                            rdc_score = rdc_score + tensor[0]
                        t2_log[rdc_score].append(hit)
                rdc_score_bins = t2_log.keys()
                rdc_score_bins.sort()
                for k in range(len(rdc_score_bins)):
                    hits = t2_log[rdc_score_bins[k]]
                    for hit in hits:
                        smotif_seq = hit[3][1]
                        if smotif_seq not in seqs:
                            seqs.append(smotif_seq)
                            reduced_dump_log.append(hit)
                            print "final sele", hit[0][1][0][0], keys[i], noe_energy_bins[j], rdc_score_bins[k]
                            count_hits += 1
                if count_hits >= num_hits:
                    break
        if count_hits >= num_hits:
            break
    if count_hits >= num_hits:
        pass
    else:
        print "could only extract ", len(reduced_dump_log), count_hits


    # io.dumpPickle(str(previous_smotif_index) + "_tophits.pickle", dump_pickle)
    io.dumpGzipPickle(str(previous_smotif_index) + "_tophits.gzip", reduced_dump_log)
    print "actual number in top hits ", len(reduced_dump_log)
    return range(len(reduced_dump_log))


def getRunSeq(num_hits, stage):
    """
    generate run seq, a seq list of pairs of
    indexes of profiles for job scheduling
    """
    map_route = []
    ss_profiles = io.readPickle("ss_profiles.pickle")
    if os.path.isfile("contacts_route.pickle"):
        map_route = io.readPickle("contacts_route.pickle")
    elif os.path.isfile("pcs_route.pickle"):
        map_route = io.readPickle("pcs_route.pickle")
    elif os.path.isfile("rdc_route.pickle"):
        map_route = io.readPickle("rdc_route.pickle")

    try:
        next_index, next_smotif = getNextSmotif(map_route)
        print next_index, next_smotif
    except TypeError:
        return [999], 999

    direction = next_smotif[-1]
    if direction == 'left':
        next_ss_list = ss_profiles[next_smotif[0]]
    else:
        next_ss_list = ss_profiles[next_smotif[1]]
    # get and make a list of top 10(n) of the previous run
    #top_hits = makeTopPickle(next_index - 1, num_hits, stage)  # send the previous Smotif index
    top_hits = makeTopPickle2(next_index - 1, num_hits, stage)  # send the previous Smotif index

    # delete two stages down pickled files
    # check_pickle = str(next_index - 2) + str("_*_*.pickle")
    check_pickle = str(next_index - 2) + str("_*_*.gzip")
    file_list = glob.glob(check_pickle)

    if len(file_list) > 10:
        remove = "rm " + check_pickle
        os.system(remove)

    if top_hits:
        run_seq = []
        for i in range(len(top_hits)):
            for j in range(len(next_ss_list)):
                run_seq.append([i, j])
        return run_seq, next_index
