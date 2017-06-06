import collections
import math

import utility.stage2_util as s2util


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
