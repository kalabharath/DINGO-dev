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
import filters.contacts.contacts_filter  as Cfilter
import filters.rmsd.qcp as qcp



def getfromDB(searched_smotifs, current_ss, direction):

    previous_ss = searched_smotifs[-1]

    #print current_ss, previous_ss

    if direction == 'left':#double check this implementation

        smotif_def = sm.getSmotif(current_ss, previous_ss)
    else:
        smotif_def = sm.getSmotif(previous_ss, current_ss)

    return sm.readSmotifDatabase(smotif_def)


def orderSSE(previous_sse, current_sse):

    previous_sse.append(current_sse)
    return previous_sse


def SmotifSearch(index_array):

    #print index_array

    preSSE = uts2.getPreviousSmotif3(index_array[0])

    preSSE = preSSE[0]

    current_ss, direction = uts2.getSS2(index_array[1])

    #print current_ss, direction
    csmotif_data = getfromDB(preSSE[-1], current_ss, direction)

    exp_data = io.readPickle("exp_data.pickle")

    dump_log = []

    sse_ordered = orderSSE(preSSE[-1], current_ss)


    for i in range(0, len(csmotif_data)):

        csmotif = csmotif_data[i]

        ##QCP RMSD

        rmsd, transformed_coos = qcp.rmsdQCP3(preSSE[0],csmotif, direction)

        if rmsd <= 2.0:

            ## Sequence filter, align native and smotif aa_seq as a measure of sequence similarity = structure similarity
            csse_seq, seq_identity, blosum62_score, bool_sequence_similarity \
            = Sfilter.S2SequenceSimilarity(current_ss, csmotif, direction, exp_data, threshold=30)
            #print csse_seq

            ## Contacts filter,
            no_of_contacts, percent_of_satisfied_contacts \
            = Cfilter.S2ContactPredicition(transformed_coos, sse_ordered, exp_data, threshold=0.8)


            if bool_sequence_similarity and percent_of_satisfied_contacts > 50.0:
                print 'blosum62 score', blosum62_score, "seq_id", seq_identity, "rmsd=", rmsd, "Contacts", percent_of_satisfied_contacts
                dump_log.append([transformed_coos,['seq_filter', csse_seq, seq_identity, blosum62_score],
                ['contacts_filter', no_of_contacts, percent_of_satisfied_contacts], sse_ordered])

    if len(dump_log) > 0 :
        io.dumpPickle("tx_"+str(index_array[0])+"_"+str(index_array[1])+".pickle",dump_log)
    return True