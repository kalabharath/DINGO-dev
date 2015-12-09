#!/usr/bin/env python

"""
Project_Name: main, File_name: stage2_search.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 7/05/15 , Time:10:05 PM

Perform stage 3x in parallel
"""
import utility.stage2_util as uts2
import utility.smotif_util as sm
import utility.io_util as io
import filters.sequence.sequence_similarity as Sfilter
import filters.pcs.pcsfilter as Pfilter
import filters.rmsd.qcp as qcp
import filters.constraints.looplengthConstraint as llc
import time


def getfromDB(previous_smotif, current_ss, direction):
    # print "previous_smotif: ", previous_smotif

    searched_smotifs = []
    for entry in previous_smotif:
        if 'smotif_def' == entry[0]:
            searched_smotifs = entry[-1]

    # ['smotif_def', [['helix', 6, 7, 5, 145, 150], ['helix', 23, 5, 1, 156, 178], ['strand', 5, 7, 8, 133, 137]]]

    if direction == 'left':
        previous_ss = searched_smotifs[0]
    else:
        previous_ss = searched_smotifs[-1]

    # print "previous_ss: ", previous_ss
    # print "current_ss : ", current_ss

    if direction == 'left':  # double check this implementation

        smotif_def = sm.getSmotif(current_ss, previous_ss)
    else:
        smotif_def = sm.getSmotif(previous_ss, current_ss)

    return sm.readSmotifDatabase(smotif_def)


def orderSSE(previous_smotif, current_sse, direction):
    """

    :param previous_smotif:
    :param current_sse:
    :return:
    """

    for entry in previous_smotif:
        # ['qcp_rmsd', transformed_coos, sse_ordered, rmsd]
        if 'qcp_rmsd' == entry[0]:
            previous_sse = entry[2]

            if direction == 'left':
                previous_sse.insert(0,current_sse)
            else:
                previous_sse.append(current_sse)
            return previous_sse


def SmotifSearch(index_array):
    # print index_array

    preSSE = uts2.getPreviousSmotif(index_array[0])
    current_ss, direction = uts2.getSS2(index_array[1])
    print current_ss, direction

    csmotif_data = getfromDB(preSSE, current_ss, direction)

    if not csmotif_data:
        # If the smotif library doesn't exist
        # Terminate further execution
        return True


    exp_data = io.readPickle("exp_data.pickle")
    exp_data_types = exp_data.keys()  # ['ss_seq', 'pcs_data', 'aa_seq', 'contacts', 'natives']
    sse_ordered = orderSSE(preSSE, current_ss, direction)
    #print sse_ordered
    dump_log = []

    no_clashes = False


    stime = time.time()
    for i in range(0, len(csmotif_data)):

        # Exclude the natives if present
        if 'natives' in exp_data_types:
            natives = exp_data['natives']
            tpdbid = csmotif_data[i][0][0]
            pdbid = tpdbid[0:4]
            if pdbid in natives:
                #Stop further execution and resume iteration
                continue

        # QCP RMSD

        rmsd, transformed_coos = qcp.rmsdQCP3(preSSE, csmotif_data[i], direction)

        if rmsd <= exp_data['rmsd']:
            loopconstraint = llc.loopConstraint(transformed_coos, sse_ordered, direction)

            if loopconstraint:
                no_clashes = qcp.clahses(transformed_coos, exp_data['clash_distance'])
            else:
                no_clashes = False

        if rmsd <= exp_data['rmsd'] and no_clashes:

            pcs_tensor_fits = []
            #print csmotif_data[i]
            tlog = []
            tlog.append(['smotif', csmotif_data[i]])
            tlog.append(['smotif_def', sse_ordered])
            tlog.append(['qcp_rmsd', transformed_coos, sse_ordered, rmsd])

            cathcodes = sm.orderCATH(preSSE, csmotif_data[i][0], direction)
            #print cathcodes
            tlog.append(['cathcodes', cathcodes])

            csse_seq, seq_identity, blosum62_score, bool_sequence_similarity \
                = Sfilter.S2SequenceSimilarity(current_ss, csmotif_data[i], direction, exp_data, threshold=40)

            concat_seq = sm.orderSeq(preSSE, csse_seq, direction)

            tlog.append(['seq_filter', concat_seq, csse_seq, seq_identity, blosum62_score])

            if 'pcs_data' in exp_data_types and seq_identity <= 30.0:
                pcs_tensor_fits = Pfilter.PCSAxRhFit2(transformed_coos, sse_ordered, exp_data, stage = 4)
                tlog.append(['PCS_filter', pcs_tensor_fits])

            if pcs_tensor_fits :
                #print csmotif_data[i][0],"seq_id", seq_identity, "rmsd=", rmsd, cathcodes
                dump_log.append(tlog)

            #Time bound search
            ctime = time.time()
            elapsed = ctime-stime
            if (elapsed/60.0)> 120.0: #stop execution after 2 hrs
                print "Breaking further execution"
                break


    if len(dump_log) > 0:
        io.dumpPickle("tx_" + str(index_array[0]) + "_" + str(index_array[1]) + ".pickle", dump_log)
    return True
