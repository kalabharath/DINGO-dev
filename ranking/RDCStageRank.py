import utility.stage2_util as s2util
import collections
import math

def rank_dump_log(dump_log, exp_data, stage):

    prob_top_hits=exp_data['prob_top_hits']
    num_hits = round(len(dump_log) * prob_top_hits[stage-1], 1)
    new_dict = collections.defaultdict(list)

    rdc_filter = False
    noe_filter = False
    for hit in dump_log:
        # thread_data contains data from each search and filter thread.
        for data_filter in hit:
            if data_filter[0] == 'RDC_filter':
                rdc_filter = True
                rdc_data = data_filter
                Nchi = s2util.rdcSumChi(rdc_data, stage)
                for filter in hit:
                    if filter[0] == 'NOE_filter':
                        noe_filter = True
                        noe_fmeasure = filter[1]
                        Nchi = Nchi / math.pow(10, noe_fmeasure)
                        new_dict[Nchi].append(hit)
                if not noe_filter:
                    new_dict[Nchi].append(hit)

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
                if ent[0] == 'smotif':
                    name = ent[1][0]
                if ent[0] == 'seq_filter':
                    seq_filter = ent
                    smotif_seq = seq_filter[1]
                if ent[0] == 'RDC_filter':
                    rdc_data = ent
                    Nchi = s2util.rdcSumChi(rdc_data, stage)
                    if noe_filter:
                        for ent in entry:
                            if ent[0] == 'NOE_filter':
                                noe_fmeasure = ent[1]
                                Nchi = Nchi / math.pow(10, noe_fmeasure)
                    else:
                        Nchi = s2util.rdcSumChi(rdc_data, stage)

            if smotif_seq not in seqs:
                seqs.append(smotif_seq)
                # non_redundant.setdefault(Nchi, []).append(entry)
                #non_redundant[Nchi].append(entry)
                reduced_dump_log.append(entry)
                count_hits += 1
        if count_hits >= num_hits:
            break
    print "Reducing the amount of data to:",prob_top_hits[stage-1], len(reduced_dump_log), len(dump_log)
    return reduced_dump_log