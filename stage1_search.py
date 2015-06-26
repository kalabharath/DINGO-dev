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
    # TODO list the complete details for how this function works

    # print index_array
    s1_def, s2_def = getSSdef(index_array)
    smotif_def = sm.getSmotif(s1_def, s2_def)
    print s1_def, s2_def

    smotif_data = sm.readSmotifDatabase(smotif_def)
    exp_data = io.readPickle("exp_data.pickle")
    exp_data_types = exp_data.keys() #['ss_seq', 'pcs_data', 'aa_seq', 'contacts']

    dump_log=[]

    for i in range(0, len(smotif_data)):

        tlog = []
        tlog.append(['smotif',smotif_data[i]])
        tlog.append(['smotif_def',[s1_def, s2_def]])


        # TODO clever use of variable names

        ### Filters

        if 'aa_seq' in exp_data_types:
            smotif_seq, seq_identity, blosum62_score  = \
                Sfilter.SequenceSimilarity(s1_def, s2_def, smotif_data[i], exp_data)
            tlog.append(['seq_filter', smotif_seq, seq_identity, blosum62_score])

        if 'contacts' in exp_data_types:
            no_of_contacts, percent_of_satisfied_contacts = \
                Cfilter.ContactPredicition(s1_def, s2_def, smotif_data[i], exp_data)
            tlog.append(['contacts_filter', no_of_contacts, percent_of_satisfied_contacts])

        if 'pcs_data' in exp_data_types:
            pcs_tensor_fits = Pfilter.PCSAxRhFit(s1_def, s2_def, smotif_data[i], exp_data)
            tlog.append(['PCS_filter', pcs_tensor_fits])

        ### Filters

        if seq_identity > 40.0 and percent_of_satisfied_contacts > 50.0 :
            # print index_array, s1_def, s2_def
            # print smotif_def, len(smotif_data)
            print smotif_data[i][0][0], 'blosum62 score', blosum62_score, \
                "seq_id", seq_identity, "i=", i, "/", len(smotif_data), percent_of_satisfied_contacts
            print pcs_tensor_fits
            dump_log.append(tlog)

    if dump_log:
        print "num of hits", len(dump_log)
        io.dumpPickle('0_'+str(index_array[0])+"_"+str(index_array[1])+".pickle",dump_log)

    return True
