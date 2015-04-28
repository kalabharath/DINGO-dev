#!/usr/bin/env python

"""
Project_Name: main/filters/pcsfilter, File_name: pcsfilter.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 24/04/15 , Time: 01:30 PM

fit PCS and apply Ax,Rh filters
"""
import utility.io_util as io

def getPCSData():
    """
    This is a very very bad implementation, do refactoring at a a later stage
    :return:
    """
    pcs_broker = io.readInputDataFiles('input_data.txt', ['pcs_broker'])
    psi_pred = io.readInputDataFiles('input_data.txt', ['psipred_file'])
    ss_seq = io.readPsiPred(psi_pred)
    pcs_data = io.getPcsTagInfo(ss_seq, pcs_broker)
    return pcs_data

def getHN(smotif, atom_type='H'):
    """
    parse the x,y,z of individual SSE into seperate arrays
    :param smotif:
    :return: arrays of coordinates of individual SSE in seperate arrays
    """
    #TODO return arrays for any given atom type(s)

    rH1 , rH2 = [], []
    for entry in  smotif[0][1]:
        # res_no, aa_type, atom_type, x, y, z
        if entry[2] == atom_type:
            x, y, z = entry[3], entry[4], entry[5]
            rH1.append([x, y, z])

    for entry in  smotif[0][2]:
        # res_no, aa_type, atom_type, x, y, z
        if entry[2] == atom_type:
            x, y, z = entry[3], entry[4], entry[5]
            rH2.append([x, y, z])
    return rH1, rH2


def match_pcss_HN():
    return  True
def PCSAxRhFit(s1_def, s2_def, smotif, threshold = 0.05):
    """

    :param s1_def:
    :param s2_def:
    :param smotif:
    :param threshold:
    :return:
    """

    ss1_list = range(s1_def[4], s1_def[5]+1)
    ss2_list = range(s2_def[4], s2_def[5]+1)

    smotif_ss1 = range(int(smotif[0][0][1]), int(smotif[0][0][2])+1 )
    smotif_ss2 = range(int(smotif[0][0][3]), int(smotif[0][0][4])+1 )

    print ss1_list, ss2_list
    print smotif_ss1, smotif_ss2
    print smotif[0][0]

    rH1, rH2 = getHN(smotif, atom_type='H')

    pcs_data = getPCSData()
    ntags = len(pcs_data)
    for tag in range(0,ntags):
        print tag

    return True