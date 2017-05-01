import collections

import utility.stage2_util as s2util


def rank_dump_log(dump_log, exp_data, stage):

    rank_top_hits=exp_data['rank_top_hits']
    num_hits = rank_top_hits[stage-1]
    new_dict = collections.defaultdict(list)
    pcs_filter = False
    rdc_filter = False
    noe_filter = False
    global_noe_filter = False
    for hit in dump_log:
        # thread_data contains data from each search and filter thread.
        for data_filter in hit:
            if data_filter[0] == 'PCS_filter':
                pcs_filter = True
                pcs_data = data_filter
                Nchi = s2util.getNchiSum(pcs_data, stage)
                # new_dict.setdefault(Nchi, []).append(entry)
                new_dict[Nchi].append(hit)
            if data_filter[0] == 'RDC_filter':
                rdc_filter = True
                rdc_data = data_filter
                Nchi = s2util.rdcSumChi(rdc_data, stage)
                new_dict[Nchi].append(hit)

    # ************************************************
    # Exclude the redundant entries and rank top hits
    # ************************************************

    keys = new_dict.keys()
    keys.sort()
    if global_noe_filter:
        keys.reverse()

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
                if ent[0] == 'smotif':
                    name = ent[1][0]
                if ent[0] == 'seq_filter':
                    seq_filter = ent
                    smotif_seq = seq_filter[1]
            if smotif_seq not in seqs:
                seqs.append(smotif_seq)
                reduced_dump_log.append(entry)
                count_hits += 1
        if count_hits >= num_hits:
            break
    print "Reducing the amount of data to:",rank_top_hits[stage-1], len(reduced_dump_log), len(dump_log)
    return reduced_dump_log