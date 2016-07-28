#!/usr/bin/env python

"""
Project_Name: main, File_name: stage2_search.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 7/05/15 , Time:10:05 PM

Perform stage 2 in parallel
"""

import filters.constraints.looplengthConstraint as llc
import filters.contacts.evfoldContacts as Evofilter
import filters.pcs.pcsfilter as Pfilter
import filters.rmsd.qcp as qcp
import filters.sequence.sequence_similarity as Sfilter
import utility.io_util as io
import utility.smotif_util as sm
import utility.stage2_util as uts2


def getfromDB(previous_smotif, current_ss, direction, database_cutoff):
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

    return sm.readSmotifDatabase(smotif_def, database_cutoff)


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
    """
    Main()
    :param index_array:
    :return:
    """

    exp_data = io.readPickle("exp_data.pickle")
    exp_data_types = exp_data.keys()  # ['ss_seq', 'pcs_data', 'aa_seq', 'contacts']

    psmotif = uts2.getPreviousSmotif(index_array[0])  # TODO refactor these to be more general in reading the route
    current_ss, direction = uts2.getSS2(index_array[1])  # TODO refactor these to be more general in reading the route

    csmotif_data = getfromDB(psmotif, current_ss, direction, exp_data['database_cutoff'])

    if not csmotif_data:
        # If the smotif library doesn't exist
        # Terminate further execution
        return True

    """
    always narrow down to previous sse and current sse and operate on them individually
    """
    sse_ordered = orderSSE(psmotif, current_ss, direction)
    dump_log = []
    no_clashes = False

    # ************************************************************************************************
    # Main
    # The 'for' loop below iterates over all of the Smotifs and applies various filters
    # This is the place to add new filters as you desire. For starters, look at Sequence filter.
    # ************************************************************************************************


    for i in range(0, len(csmotif_data)):

        # ************************************************
        # Applying different filters for the Smotif assembly
        # ************************************************

        # Exclude natives if needed
        if 'natives' in exp_data_types:
            natives = exp_data['natives']
            tpdbid = csmotif_data[i][0][0]
            pdbid = tpdbid[0:4]
            if pdbid in natives:
                # Stop further execution and
                continue

        # ************************************************
        # RMSD filter using QCP method
        # quickly filters non-overlapping smotifs
        # ************************************************

        rmsd, transformed_coos = qcp.rmsdQCP(psmotif[0], csmotif_data[i], direction)

        if rmsd <= exp_data['rmsd_cutoff'][1]:
            # Loop constraint restricts the overlapping smotifs is not drifted far away.
            loop_constraint = llc.loopConstraint(transformed_coos, sse_ordered, direction)
            if loop_constraint:
                # Check whether the SSEs with in the assembled smotifs are clashing to one another
                no_clashes = qcp.clahses(transformed_coos, exp_data['clash_distance'])
            else:
                no_clashes = False


        if rmsd <= exp_data['rmsd_cutoff'][1] and no_clashes:
            # Prepare temp log array to save data at the end
            tlog = []
            pcs_tensor_fits = []
            contact_fmeasure = []

            tlog.append(['smotif', csmotif_data[i]])
            tlog.append(['smotif_def', sse_ordered])

            tlog.append(['qcp_rmsd', transformed_coos, sse_ordered, rmsd])
            cathcodes = sm.orderCATH(psmotif, csmotif_data[i][0], direction)
            tlog.append(['cathcodes', cathcodes])

            # ************************************************
            # Sequence filter
            # Aligns the smotif seq to target seq and calculates
            # sequence identity and the alignment score
            # ************************************************

            csse_seq, seq_identity, blosum62_score = Sfilter.S2SequenceSimilarity(current_ss, csmotif_data[i],
                                                                                  direction, exp_data)

            # concat current to previous seq
            concat_seq = sm.orderSeq(psmotif, csse_seq, direction)

            tlog.append(['seq_filter', concat_seq, csse_seq, seq_identity, blosum62_score])

            # ************************************************
            # Pseudocontact Shift filter
            # uses experimental PCS data to filter Smotifs
            # scoring based on normalised chisqr
            # ************************************************

            if 'pcs_data' in exp_data_types and seq_identity >= 10.0:
                pcs_tensor_fits = Pfilter.PCSAxRhFit2(transformed_coos, sse_ordered, exp_data, stage=2)
                tlog.append(['PCS_filter', pcs_tensor_fits])

            # ************************************************
            # Contacts filter
            # uses the contact data obtained from EVfold server
            # tp score a given smotif
            # ************************************************
            if 'contact_matrix' in exp_data_types:

                contact_fmeasure, plm_score = Evofilter.s2EVcouplings(transformed_coos, sse_ordered,
                                                                      exp_data['contact_matrix'],
                                                                      exp_data['plm_scores'],
                                                                      contacts_cutoff=9.0)
                if contact_fmeasure and plm_score:

                    if contact_fmeasure >= 0.6:
                        contact_score = (contact_fmeasure * 2) + (plm_score * 0.1) + (seq_identity * (0.01) * (5))
                        # print csmotif_data[i][
                        #    0], 'fmeasure', contact_fmeasure, "seq_id", seq_identity, "rmsd=", rmsd, cathcodes
                    elif contact_fmeasure > 0.5 and contact_fmeasure < 0.6:
                        contact_score = (contact_fmeasure) + (plm_score * 0.1) + (seq_identity * (0.01) * (5))
                        # continue
                    else:
                        contact_score = (contact_fmeasure) + (plm_score * 0.1) + (seq_identity * (0.01) * (5))
                        # continue
                    tlog.append(['Evofilter', contact_score])

            if pcs_tensor_fits or contact_fmeasure:
                #dump data to the disk
                # print csmotif_data[i][0], 'blosum62 score', blosum62_score, "seq_id", seq_identity, "rmsd=", rmsd, cathcodes
                dump_log.append(tlog)

    # prevent dumping empty arrays with no data
    if len(dump_log) > 0:
        print "num of hits", len(dump_log),
        io.dumpPickle("tx_" + str(index_array[0]) + "_" + str(index_array[1]) + ".pickle", dump_log)

    return True
