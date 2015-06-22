#!/usr/bin/env python

"""
Project_Name: main/utility, File_name: smotif_util.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 14/04/15 , Time:12:54 PM
"""

import utility.io_util as io


def getSmotif(s1, s2):
    """
    Return Smotif based on SS definition

    ('helix-helix', 6, 23, 4)
    new_smotif_file = hh_6_23.db
    ['helix', 6, 9, 4, 146, 151] ['helix', 23, 4, 1, 156, 178]
    """
    if s1[0] == 'helix':
        s1_type = 'h'
    if s1[0] == 'strand':
        s1_type = 's'
    if s2[0] == 'helix':
        s2_type = 'h'
    if s2[0] == 'strand':
        s2_type = 's'
    smotif_type = s1_type + s2_type
    s1_len = s1[1]
    s2_len = s2[1]
    smotif = [smotif_type, s1_len, s2_len]
    return smotif


def readSmotifDatabase(smotif):
    """

    :param smotif:
    :return:
    """
    # TODO option to parse in database path

    #smotif_db_path = "/home/kalabharath/zinr/smotif_database/"
    smotif_db_path = "/home/kalabharath/zinr/main/smotif_cen_db/"
    file_name = smotif[0] + "_" + str(smotif[1]) + "_" + str(smotif[2]) + ".db"
    fin = smotif_db_path + file_name
    smotif_data = io.readPickle(fin)

    return smotif_data


