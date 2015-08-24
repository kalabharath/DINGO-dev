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
import filters.contacts.contacts_filter  as Cfilter
import filters.pcs.pcsfilter as Pfilter
import filters.rmsd.qcp as qcp


def getfromDB(previous_smotif, current_ss, direction):

    for entry in previous_smotif:
        if 'smotif_def' == entry[0]:
            psmotif = entry[-1]

    if direction == 'left':
        previous_ss = psmotif[0]
    else:
        previous_ss = psmotif[1]

    # current_ss, previous_ss
    if direction == 'left':#double check this implementation
        smotif_def = sm.getSmotif(current_ss, previous_ss)
    else:
        smotif_def = sm.getSmotif(previous_ss, current_ss)

    return sm.readSmotifDatabase(smotif_def)


def orderSSE(previous_smotif, current_sse):

    previous_seq = []
    for entry in previous_smotif:
        if 'smotif_def' == entry[0]:
            previous_seq = entry

    ordered_SSE = []
    for sse_array in previous_seq[1]:
        ordered_SSE.append(sse_array)
    ordered_SSE.append(current_sse)
    return ordered_SSE


def SmotifSearch(index_array):

    #print index_array

    psmotif = uts2.getPreviousSmotif(index_array[0])

    current_ss, direction = uts2.getSS2(index_array[1])
    csmotif_data = getfromDB(psmotif, current_ss, direction)

    exp_data = io.readPickle("exp_data.pickle")
    exp_data_types = exp_data.keys()  # ['ss_seq', 'pcs_data', 'aa_seq', 'contacts']

    """
    always narrow down to previous sse and current sse and operate on them individually

    """
    sse_ordered = orderSSE(psmotif, current_ss)

    dump_log = []

    for i in range(0, len(csmotif_data)):


        # QCP RMSD
        rmsd, transformed_coos = qcp.rmsdQCP(psmotif[0],csmotif_data[i], direction)
        no_clashes = qcp.clahses(transformed_coos)

        #if rmsd <= 1.5 and no_clashes:
        if csmotif_data[i][0][0] == '2z2iA00':
            tlog = []
            tlog.append(['smotif', csmotif_data[i]])
            tlog.append(['smotif_def', sse_ordered])
            tlog.append(['qcp_rmsd', transformed_coos, sse_ordered, rmsd])

            ## Sequence filter, align native and smotif aa_seq as a measure of sequence similarity = structure similarity

            if 'aa_seq' in exp_data_types:
                csse_seq, seq_identity, blosum62_score, bool_sequence_similarity \
                    = Sfilter.S2SequenceSimilarity(current_ss, csmotif_data[i], direction, exp_data, threshold=40)
                tlog.append(['seq_filter', csse_seq, seq_identity, blosum62_score])

            if 'contacts' in exp_data_types:
                no_of_contacts, percent_of_satisfied_contacts \
                    = Cfilter.S2ContactPredicition(transformed_coos, sse_ordered, exp_data)
                tlog.append(['contacts_filter', no_of_contacts, percent_of_satisfied_contacts])

            if 'pcs_data' in exp_data_types:
                pcs_tensor_fits = Pfilter.PCSAxRhFit2(transformed_coos, sse_ordered, exp_data)
                tlog.append(['PCS_filter', pcs_tensor_fits])

            #if pcs_tensor_fits and seq_identity > 40:
            if True:
                print "rmsd", rmsd
                print csmotif_data[i][0]
                print pcs_tensor_fits
                print 'blosum62 score', blosum62_score, "seq_id", seq_identity, "rmsd=", rmsd, "Contacts", percent_of_satisfied_contacts
                dump_log.append(tlog)

    if len(dump_log) > 0 :
        io.dumpPickle("tx_"+str(index_array[0])+"_"+str(index_array[1])+".pickle",dump_log)

    return True
