import collections
import math

import utility.stage2_util as s2util

def limitParents(current_parent, smotif_parents):
    count = 0
    if len(smotif_parents) == 0:
        return True
    for sp in smotif_parents:
        if sp == current_parent:
            count += 1
            if count > 5:
                return False
    return True

def rank_assembly(dump_log, num_hits):


    new_dict = collections.defaultdict(list)

    """
    0 smotif
    1 smotif_def
    2 cathcodes
    3 qcp_rmsd
    4 seq_filter
    5 NOE_filter
    6 RDC_filter
    6 Ref_RMSD
    """

    for hit in dump_log:
        # thread_data contains data from each search and filter thread.
        # initialize total score array
        #if hit[4][0] == 'NOE_filter':
        noe_energy = hit[5][3]
        noe_energy = round(noe_energy, 3)
        new_dict[noe_energy].append(hit)

    keys = new_dict.keys()
    keys.sort()
    # Rank based on NOE energy

    reduced_dump_log = []
    seqs = []
    parents = []
    count_hits = 0

    for i in range(len(keys)):
        entries = new_dict[keys[i]]

        if count_hits >= num_hits:
            break

        if len(entries) == 1:
            smotif_seq = entries[0][4][1]
            smotif_parents = entries[0][3][2]
            if (smotif_seq not in seqs) and (limitParents(smotif_parents, parents)):
                seqs.append(smotif_seq)
                parents.append(smotif_parents)
                print parents
                reduced_dump_log.append(entries[0])
                count_hits += 1
                print "final sele", entries[0][0][1][0][0], keys[i]
                if count_hits >= num_hits:
                    break
        else:
            t2_log = collections.defaultdict(list)
            for hit in entries:
                #if hit[5][0] == 'RDC_filter':
                rdc_tensors = hit[6][1]
                rdc_score = 0
                for tensor in rdc_tensors:
                    rdc_score = rdc_score + tensor[0]
                t2_log[rdc_score].append(hit)
            rdc_score_bins = t2_log.keys()
            rdc_score_bins.sort()
            for k in range(len(rdc_score_bins)):
                hits = t2_log[rdc_score_bins[k]]
                for hit in hits:
                    smotif_seq = hit[4][1]
                    smotif_parents = hit[3][2]
                    #if (smotif_seq not in seqs) and (smotif_parents not in parents):
                    if (smotif_seq not in seqs) and (limitParents(smotif_parents, parents)):
                        seqs.append(smotif_seq)
                        parents.append(smotif_parents)
                        print parents
                        reduced_dump_log.append(hit)
                        print "final sele", hit[0][1][0][0], keys[i], rdc_score_bins[k]
                        count_hits += 1
                    if count_hits >= num_hits:
                        break
                if count_hits >= num_hits:
                    break
            if count_hits >= num_hits:
                break
    if count_hits >= num_hits:
        pass
    else:
        print "could only extract ", len(reduced_dump_log), count_hits

    return reduced_dump_log

def rank_assembly_old(dump_log, exp_data, stage):
    rank_top_hits = exp_data['rank_top_hits']
    num_hits = rank_top_hits[stage - 1]

    new_dict = collections.defaultdict(list)

    """
    0 smotif
    1 smotif_def
    2 cathcodes
    3 seq_filter
    4 NOE_filter
    5 RDC_filter
    6 Ref_RMSD
    """

    for hit in dump_log:
        # thread_data contains data from each search and filter thread.
        # initialize total score array
        #if hit[4][0] == 'NOE_filter':
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
        if len(
                entries) == 1:  # There is only one entry in this no_of_noes bin just check of existing sequences and move on
            smotif_seq = entries[0][3][1]
            if smotif_seq not in seqs:
                seqs.append(smotif_seq)
                reduced_dump_log.append(entries[0])
                count_hits += 1
        else:
            t_log = collections.defaultdict(list)
            for hit in entries:  # filter on noe_energy
                #if hit[4][0] == 'NOE_filter':
                noe_energy = hit[4][3]
                noe_energy = round(noe_energy, 2)
                t_log[noe_energy].append(hit)
            noe_energy_bins = t_log.keys()
            noe_energy_bins.sort()

            for j in range(len(noe_energy_bins)):  # filter on RDC score
                t2_log = collections.defaultdict(list)
                hits = t_log[noe_energy_bins[j]]
                for hit in hits:
                    #if hit[5][0] == 'RDC_filter':
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
                            count_hits += 1
                if count_hits >= num_hits:
                    break
        if count_hits >= num_hits:
            break
        else:
            pass
    print "Reducing the amount of data to:", rank_top_hits[stage - 1], len(reduced_dump_log), len(dump_log)
    return reduced_dump_log



def rank_dump_log(dump_log, exp_data, stage):
    rank_top_hits = exp_data['rank_top_hits']
    num_hits = rank_top_hits[stage - 1]
    new_dict = collections.defaultdict(list)
    rdc_constant = 0.0
    for hit in dump_log:
        # thread_data contains data from each search and filter thread.
        # initialize total score array
        total_score = {}
        for data_filter in hit:

            if data_filter[0] == 'PCS_filter':
                pcs_data = data_filter
                pcsscore = s2util.getNchiSum(pcs_data, stage)
                total_score['pcs_score'] = pcsscore

            if data_filter[0] == 'RDC_filter':
                rdc_data = data_filter
                #Nchi = s2util.rdcSumChi(rdc_data, stage)
                log_likelihood = data_filter[2]
                rdc_tensors = data_filter[1]
                for tensor in rdc_tensors:
                    rdc_constant = rdc_constant + tensor[0]
                rdc_constant = rdc_constant * 1e-10
                total_score['rdc_score'] = log_likelihood

            if data_filter[0] == 'NOE_filter':
                noe_probability = data_filter[1]
                log_likelihood = -(math.log(noe_probability))
                total_score['noe_score'] = log_likelihood

                # calculate the total score and append the hit
        if total_score:
            keys = total_score.keys()
            keys = ['noe_score','rdc_score']
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
    non_redundant = collections.defaultdict(list)
    reduced_dump_log = []
    seqs = []
    smotif_seq = ''
    Nchi = 0.0
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
                reduced_dump_log.append(entry)
                count_hits += 1
                if count_hits >= num_hits:
                    break
            if count_hits >= num_hits:
                break
        if count_hits >= num_hits:
            break
    print "Reducing the amount of data to:", rank_top_hits[stage - 1], len(reduced_dump_log), len(dump_log)
    return reduced_dump_log
