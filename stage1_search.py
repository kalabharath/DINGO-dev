#!/usr/bin/env python

"""
Project_Name: main, File_name: stage_1_search.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 14/04/15 , Time:01:05 PM

stage 1 in parallel
"""

import utility.stage1_util as uts1
import utility.smotif_util as sm
import utility.io_util as io
import filters.sequence.sequence_similarity as Sfilter
import filters.pcs.pcsfilter as Pfilter
import time


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
    ctime = time.time()

    for i in range(0, len(smotif_data)):

        stime = time.time()
        elapsed = ctime-stime
        if elapsed/60.0> 120.0: #stop execution after 2 hrs
            return True


        # Save CPU time by excluding natives
        if 'natives' in exp_data_types:
            natives = exp_data['natives']
            tpdbid = smotif_data[i][0][0]
            pdbid = tpdbid[0:4]
            if pdbid in natives:
                #Stop further execution and iterate
                continue
        tlog = []
        pcs_tensor_fits = []
        tlog.append(['smotif',smotif_data[i]])
        tlog.append(['smotif_def',[s1_def, s2_def]])
        tlog.append(['cathcodes',[smotif_data[i][0]]])

        ### Filters
        # TODO clever use of variable names
        # Add new modules here

        #Calculate Sequence identity first

        smotif_seq, seq_identity, blosum62_score  = \
                Sfilter.SequenceSimilarity(s1_def, s2_def, smotif_data[i], exp_data)
        tlog.append(['seq_filter', smotif_seq, seq_identity, blosum62_score])

        if 'pcs_data' in exp_data_types and seq_identity >= 0.0:
            pcs_tensor_fits = Pfilter.PCSAxRhFit(s1_def, s2_def, smotif_data[i], exp_data)
            tlog.append(['PCS_filter', pcs_tensor_fits])


        if pcs_tensor_fits:
        #if blosum62_score > 0.0 :
            # print index_array, s1_def, s2_def
            # print smotif_def, len(smotif_data)
            print smotif_data[i][0][0], 'blosum62 score', blosum62_score, \
                "seq_id", seq_identity, "i=", i, "/", len(smotif_data)
            #print pcs_tensor_fits
            dump_log.append(tlog)
                #Time bound search

        stime = time.time()
        elapsed = ctime-stime
        if elapsed/60.0> 120.0: #stop execution after 2 hrs
            if len(dump_log) > 0:
                io.dumpPickle("tx_" + str(index_array[0]) + "_" + str(index_array[1]) + ".pickle", dump_log)
                return True
            else:
                return True
    if dump_log:
        print "num of hits", len(dump_log)
        io.dumpPickle('0_'+str(index_array[0])+"_"+str(index_array[1])+".pickle",dump_log)

    return True
