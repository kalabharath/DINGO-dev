#!/usr/bin/env python

"""
Project_Name: main, File_name: stage_1_search.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 14/04/15 , Time:01:05 PM

Perform stage 1 in perfect parallel
"""

import utility.stage1_util as uts1
import utility.smotif_util as sm
#import filters.contacts.


def getSSdef(index_array):
    """

    :param index_array:
    :return:
    """
    s1_list, s2_list = uts1.getSSlist()
    return s1_list[index_array[0]], s2_list[index_array[1]]



def SmotifSearch(index_array):
    """
    #TODO list the complete details for how this function works
    :param index_array:
    :return:
    """

    s1_def, s2_def = getSSdef(index_array)
    print index_array
    print s1_def
    print s2_def
    smotif = sm.getSmotif(s1_def, s2_def)
    print smotif
    smotif_data = sm.readSmotifDatabase(smotif)
    print len(smotif_data)

    #for smotif in smotif_data:
    for i in range(0,10):
        print smotif_data[i][0][0]
        smotif = smotif_data[i]
        print smotif

        # Apply various filters. Nested filters may be a bad idea
        # TODO explore the idea of using nested filters
        sequence_similarity = filter.SequenceSimilarity(s1_def, s2_def, smotif, threshold = 0.3)
        contacts_predicition = filter.ContactPredicition(s1_def, s2_def, smotif,threshold = 0.7)
        pcs_axrh_fit_filters = filter.PCSAxRhFit(s1_def, s2_def, threshold = 0.05)








    return  True

