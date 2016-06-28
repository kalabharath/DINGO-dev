#!/usr/bin/env python

"""
Project_Name: main, File_name: stage_1_search.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 14/04/15 , Time:01:05 PM

stage 1 in parallel
"""

import filters.contacts.evfoldContacts as Evofilter
import filters.pcs.pcsfilter as Pfilter
import filters.sequence.sequence_similarity as Sfilter
import utility.io_util as io
import utility.smotif_util as sm
import utility.stage1_util as uts1


def getSSdef(index_array):
    """

    :param index_array:
    :return:
    """
    s1_list, s2_list = uts1.getSSlist()
    return s1_list[index_array[0]], s2_list[index_array[1]]


def SmotifSearch(index_array):
    """
    Main ()
    :param index_array:
    :return:
    """

    # print index_array
    s1_def, s2_def = getSSdef(index_array)
    smotif_def = sm.getSmotif(s1_def, s2_def)
    # print s1_def, s2_def

    exp_data = io.readPickle("exp_data.pickle")
    exp_data_types = exp_data.keys()  # ['ss_seq', 'pcs_data', 'aa_seq', 'contacts']

    smotif_data = sm.readSmotifDatabase(smotif_def, exp_data['database_cutoff'])

    if not smotif_data:
        # If the smotif library doesn't exist, terminate further execution.
        return True

    dump_log = []
    contact_fmeasure = []

    for i in range(0, len(smotif_data)):
        # loop over for all of the entries in the smotif_db file

        # ************************************************
        # Excluding the natives
        # ************************************************

        if 'natives' in exp_data_types:
            natives = exp_data['natives']
            tpdbid = smotif_data[i][0][0]
            pdbid = tpdbid[0:4]
            if pdbid in natives:
                # Stop further execution, but, iterate.
                continue

        # ************************************************
        # Applying different filters to Smotifs
        # Prepare temp log array to save data at the end
        # ************************************************

        tlog = []
        pcs_tensor_fits = []
        tlog.append(['smotif', smotif_data[i]])
        tlog.append(['smotif_def', [s1_def, s2_def]])
        tlog.append(['cathcodes', [smotif_data[i][0]]])

        # ************************************************
        # Sequence filter
        # Aligns the smotif seq to target seq and calculates
        # sequence identity and the alignment score
        # ************************************************

        smotif_seq, seq_identity, blosum62_score = \
            Sfilter.SequenceSimilarity(s1_def, s2_def, smotif_data[i], exp_data)
        tlog.append(['seq_filter', smotif_seq, seq_identity, blosum62_score])

        # ************************************************
        # Contacts filter
        # uses the contact data obtained from EVfold server
        # tp score a given smotif
        # ************************************************


        if 'contact_matrix' in exp_data_types:
            contact_fmeasure, plm_score = Evofilter.s1EVcouplings(s1_def, s2_def, smotif_data[i],
                                                                  exp_data['contact_matrix'],
                                                                  exp_data['plm_scores'],
                                                                  contacts_cutoff=7.0)
            print contact_fmeasure
            if contact_fmeasure and plm_score:

                if contact_fmeasure >= 0.6:
                    contact_score = (contact_fmeasure * 2) + (plm_score * 0.1) + (seq_identity * (0.01) * (5))
                elif contact_fmeasure > 0.5 and contact_fmeasure < 0.6:
                    contact_score = contact_fmeasure + (plm_score * 0.1) + (seq_identity * (0.01) * (5))
                else:
                    continue
                tlog.append(['Evofilter', contact_score])

        # ************************************************
        # Pseudocontact Shift filter
        # uses experimental PCS data to filter Smotifs
        # scoring based on normalised chisqr
        # ************************************************

        if 'pcs_data' in exp_data_types and seq_identity >= 0.0:
            pcs_tensor_fits = Pfilter.PCSAxRhFit(s1_def, s2_def, smotif_data[i], exp_data)
            tlog.append(['PCS_filter', pcs_tensor_fits])

        # Dump the data to the disk
        if pcs_tensor_fits or contact_fmeasure:
            # print smotif_data[i][0][0], "seq_id", seq_identity, "i=", i, "/", len(smotif_data)
            dump_log.append(tlog)

    # Save all of the hits in pickled arrays
    if dump_log:
        print "num of hits", len(dump_log)
        io.dumpPickle('0_' + str(index_array[0]) + "_" + str(index_array[1]) + ".pickle", dump_log)

    return True
