#!/usr/bin/env python

"""
Project_Name: main, File_name: stage_1_search.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 14/04/15 , Time:01:05 PM

Perform stage 1 in perfect parallel
"""

import utility.stage1_util as uts1
import utility.smotif_util as sm
import filters.sequence.sequence_similarity as Sfilter
import filters.contacts.contacts_filter  as Cfilter
import utility.io_util as io
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
    #TODO list the complete details for how this function works
    :param index_array:
    :return:
    """

    s1_def, s2_def = getSSdef(index_array)
    smotif_def = sm.getSmotif(s1_def, s2_def)
    smotif_data = sm.readSmotifDatabase(smotif_def)

    exp_data = io.readPickle("exp_data.pickle")

    for i in range(0, len(smotif_data)):
    #for i in range(0,100):
        #print smotif_data[i][0][0]
        smotif = smotif_data[i]

        # TODO explore the idea of using nested filters
        # TODO clever use of variable names, could be confusing if someone other than me works on it or confusing to me itself

        seq_id, seq_similar_score, bool_sequence_similarity = Sfilter.SequenceSimilarity(s1_def, s2_def, smotif, exp_data, threshold=41)

        contacts_predicition = Cfilter.ContactPredicition(s1_def, s2_def, smotif, exp_data, threshold=0.8)

        tensor_fits = Pfilter.PCSAxRhFit(s1_def, s2_def, smotif, exp_data, threshold=0.05)


        if bool_sequence_similarity and contacts_predicition > 50.0 :
            print index_array, s1_def, s2_def
            print smotif_def, len(smotif_data)
            print smotif[0][0], 'score', seq_similar_score, "seq_id", seq_id, "i=", i, "/", len(smotif_data), contacts_predicition
            print tensor_fits
    return True

