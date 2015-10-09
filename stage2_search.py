#!/usr/bin/env python

"""
Project_Name: main, File_name: stage2_search.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 7/05/15 , Time:10:05 PM

Perform stage 2 in parallel
"""
import utility.stage2_util as uts2
import utility.smotif_util as sm
import utility.io_util as io
import filters.sequence.sequence_similarity as Sfilter
import filters.pcs.pcsfilter as Pfilter
import filters.rmsd.qcp as qcp
import time


def getfromDB(previous_smotif, current_ss, direction):
    """
    :param previous_smotif:
    :param current_ss:
    :param direction:
    :return:
    """

    for entry in previous_smotif:
        if 'smotif_def' == entry[0]:
            psmotif = entry[-1]

    if direction == 'left':
        previous_ss = psmotif[0]
    else:
        previous_ss = psmotif[1]

    # current_ss, previous_ss
    if direction == 'left':  # double check this implementation
        smotif_def = sm.getSmotif(current_ss, previous_ss)
    else:
        smotif_def = sm.getSmotif(previous_ss, current_ss)

    return sm.readSmotifDatabase(smotif_def)


def orderSSE(previous_smotif, current_sse, direction):
    """
    :param previous_smotif:
    :param current_sse:
    :param direction:
    :return:
    """

    previous_seq = []
    for entry in previous_smotif:
        if 'smotif_def' == entry[0]:
            previous_seq = entry

    ordered_SSE = []
    for sse_array in previous_seq[1]:
        ordered_SSE.append(sse_array)

    if direction == 'left':
        ordered_SSE.insert(0, current_sse)
    else:
        ordered_SSE.append(current_sse)
    return ordered_SSE


def SmotifSearch(index_array):
    # print index_array

    psmotif = uts2.getPreviousSmotif(index_array[0])

    current_ss, direction = uts2.getSS2(index_array[1])
    csmotif_data = getfromDB(psmotif, current_ss, direction)

    if not csmotif_data:
        # If the smotif library doesn't exist
        # Terminate further execution
        return True


    exp_data = io.readPickle("exp_data.pickle")
    exp_data_types = exp_data.keys()  # ['ss_seq', 'pcs_data', 'aa_seq', 'contacts']

    """
    always narrow down to previous sse and current sse and operate on them individually

    """
    sse_ordered = orderSSE(psmotif, current_ss, direction)

    dump_log = []

    stime = time.time()

    for i in range(0, len(csmotif_data)):

        # Exclude natives if needed
        if 'natives' in exp_data_types:
            natives = exp_data['natives']
            tpdbid = csmotif_data[i][0][0]
            pdbid = tpdbid[0:4]
            if pdbid in natives:
                # Stop further execution and
                continue

        # QCP RMSD
        rmsd, transformed_coos = qcp.rmsdQCP(psmotif[0], csmotif_data[i], direction)
        #no_clashes = qcp.clahses(transformed_coos)
        no_clashes = qcp.clahses(transformed_coos, exp_data['clash_distance'])


        if rmsd <= exp_data['rmsd'] and no_clashes:
            tlog = []
            pcs_tensor_fits = []

            tlog.append(['smotif', csmotif_data[i]])
            tlog.append(['smotif_def', sse_ordered])
            tlog.append(['qcp_rmsd', transformed_coos, sse_ordered, rmsd])
            cathcodes = sm.orderCATH(psmotif, csmotif_data[i][0], direction)


            tlog.append(['cathcodes', cathcodes])

            ## Sequence filter, align native and smotif aa_seq as a measure of sequence similarity = structure similarity

            csse_seq, seq_identity, blosum62_score, bool_sequence_similarity \
                = Sfilter.S2SequenceSimilarity(current_ss, csmotif_data[i], direction, exp_data, threshold=40)
            # concat current to previous seq
            concat_seq = sm.orderSeq(psmotif, csse_seq, direction)

            tlog.append(['seq_filter', concat_seq, csse_seq, seq_identity, blosum62_score])

            if 'pcs_data' in exp_data_types and seq_identity >= 0.0:
                pcs_tensor_fits = Pfilter.PCSAxRhFit2(transformed_coos, sse_ordered, exp_data)
                tlog.append(['PCS_filter', pcs_tensor_fits])

            if pcs_tensor_fits:
                #print csmotif_data[i][0], 'blosum62 score', blosum62_score, "seq_id", seq_identity, "rmsd=", rmsd, cathcodes
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
