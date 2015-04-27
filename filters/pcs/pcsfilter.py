#!/usr/bin/env python

"""
Project_Name: main/filters/pcsfilter, File_name: pcsfilter.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 24/04/15 , Time: 01:30 PM

fit PCS and apply Ax,Rh filters
"""



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
    return True