#!/usr/bin/env python

"""
Project_Name: main, File_name: stage2_search.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 7/05/15 , Time:10:05 PM

Perform stage 2 in perfect parallel
"""
import utility.stage2_util as uts2
import utility.smotif_util as sm
import utility.io_util as io
import filters.sequence.sequence_similarity as Sfilter
import filters.contacts.contacts_filter  as Cfilter
import filters.rmsd.pymol_rmsd as pyrmsd
import filters.pcs.pcsfilter as Pfilter


def getfromDB(searched_smotifs, current_ss, direction):

    psmotif = searched_smotifs[-1]

    if direction == 'left':
        previous_ss = psmotif[0]
    else:
        previous_ss = psmotif[1]

    if direction == 'left':#double check this implementation
        smotif_def = sm.getSmotif(current_ss, previous_ss)
    else:
        smotif_def = sm.getSmotif(previous_ss, current_ss)

    return sm.readSmotifDatabase(smotif_def)


def SmotifSearch(index_array):

    print index_array

    psmotif = uts2.getPreviousSmotif(index_array[0])
    current_ss, direction = uts2.getSS2(index_array[1])


    csmotif_data = getfromDB(psmotif[-1], current_ss, direction)
    exp_data = io.readPickle("exp_data.pickle")
    """
    always narrow down to previous sse and current sse and operate on them individually

    """
    dump_log = []
    for i in range(0, len(csmotif_data)):
    #for i in range(0,100):
        csmotif = csmotif_data[i]
        ##Pymol RMSD

        #rmsd, transformed_coos = pyrmsd.pymol_filter(psmotif,csmotif,direction)
        bool_temp = pyrmsd.pymol_filter(psmotif[0],csmotif,direction)

        ##

        ## Sequence filter, align native and smotif aa_seq as a measure of sequence similarity = structure similarity
        smotif_seq, seq_identity, blosum62_score, bool_sequence_similarity \
            = Sfilter.S2SequenceSimilarity(current_ss, csmotif, direction, exp_data, threshold=40)
        ##

        ## Contacts filter,
        #no_of_contacts, percent_of_satisfied_contacts \
        #    = Cfilter.S2ContactPredicition(s1_def, s2_def, smotif, exp_data, threshold=0.8)

        ##
        #pcs_tensor_fits = Pfilter.PCSAxRhFit(s1_def, s2_def, smotif, exp_data, threshold=0.05)

    """
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
    """
    return True