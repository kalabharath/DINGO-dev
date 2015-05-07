#!/usr/bin/env python

"""
Project_Name: main, File_name: stage_1_search.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 14/04/15 , Time:01:05 PM

Perform stage 1 in perfect parallel
"""

import utility.stage1_util as uts1
import utility.smotif_util as sm
import utility.io_util as io
import filters.sequence.sequence_similarity as Sfilter
import filters.contacts.contacts_filter  as Cfilter
import filters.pcs.pcsfilter as Pfilter


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
    #TODO list the complete details for how this function works
    #print index_array
    s1_def, s2_def = getSSdef(index_array)
    smotif_def = sm.getSmotif(s1_def, s2_def)
    smotif_data = sm.readSmotifDatabase(smotif_def)

    exp_data = io.readPickle("exp_data.pickle")
    """
    dump_log = [[['smotif'],['seq_filter', 'smotif_seq', 'seq_identity', "blosum62_score"],
                ['contacts_filter','no_of_contacts', '%_of_contacts_observed'],
                ['PCS_filter', 'tensor_fits']]]
    """
    dump_log = []
    for i in range(0, len(smotif_data)):
    #for i in range(0,100):
        #print smotif_data[i][0][0]
        smotif = smotif_data[i]
        #print smotif[0][0]

        # TODO explore the idea of using nested filters, may not be useful!
        # TODO clever use of variable names

        smotif_seq, seq_identity, blosum62_score, bool_sequence_similarity = Sfilter.SequenceSimilarity(s1_def, s2_def, smotif, exp_data, threshold=41)
        no_of_contacts, percent_of_satisfied_contacts = Cfilter.ContactPredicition(s1_def, s2_def, smotif, exp_data, threshold=0.8)
        pcs_tensor_fits = Pfilter.PCSAxRhFit(s1_def, s2_def, smotif, exp_data, threshold=0.05)


        if bool_sequence_similarity and percent_of_satisfied_contacts > 50.0 :
            #print index_array, s1_def, s2_def
            #print smotif_def, len(smotif_data)
            #print smotif[0][0], 'blosum62 score', blosum62_score, "seq_id", seq_identity, "i=", i, "/", len(smotif_data), percent_of_satisfied_contacts
            #print pcs_tensor_fits
            dump_log.append([smotif,['seq_filter', smotif_seq, seq_identity, blosum62_score],
                             ['contacts_filter', no_of_contacts, percent_of_satisfied_contacts],
                             ['PCS_filter', pcs_tensor_fits]])
    if len(dump_log) > 1 :
        io.dumpPickle('0_'+str(index_array[0])+"_"+str(index_array[1])+".pickle",dump_log)
    return True
