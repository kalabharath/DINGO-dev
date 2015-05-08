#!/usr/bin/env python

"""
Project_Name: main/filters/sequence, File_name: sequence_similarity.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 16/04/15 , Time: 04:52 PM

Whether the Smotifs satisfy the input predicted contacts
"""

def get_distance(coo1, coo2):
    import math
    x1,y1,z1=coo1[0],coo1[1],coo1[2]
    x2,y2,z2=coo2[0],coo2[1],coo2[2]
    return math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))


def ContactPredicition(s1_def, s2_def, smotif, exp_data, threshold):

    """

    :param s1_def:
    :param s2_def:
    :param smotif:
    :param threshold:
    :return:
    """
    #     print s1_def, s2_def ['strand', 9, 4, 4, 89, 97] ['strand', 9, 4, 4, 102, 110]

    ss1_list = range(s1_def[4], s1_def[5]+1)
    ss2_list = range(s2_def[4], s2_def[5]+1)

    smotif_ss1 = range(int(smotif[0][0][1]), int(smotif[0][0][2])+1 )
    smotif_ss2 = range(int(smotif[0][0][3]), int(smotif[0][0][4])+1 )

    # print ss1_list, ss2_list
    # print smotif_ss1, smotif_ss2
    # print smotif[0][0]

    contacts = exp_data['contacts']

    residue_list = []
    contacts_in_smotif = 0
    for contact in contacts:
        if contact[0] in ss1_list and contact[1] in ss2_list:
            contacts_in_smotif +=1
            res1_index = ss1_list.index(contact[0])
            res2_index = ss2_list.index(contact[1])
            #print contact, res1_index, res2_index

            smotif_res1 = smotif_ss1[res1_index]
            smotif_res2 = smotif_ss2[res2_index]
            #print smotif_res1, smotif_res2

            for res in smotif[0][1]:
                if res[2] == 'CA' and res[0] == smotif_res1 :
                   coo1 = [res[3], res[4], res[5]]

            for res in smotif[0][2]:
                if res[2] == 'CA' and res[0] == smotif_res2 :
                   coo2 = [res[3], res[4], res[5]]
            if coo1 and coo2 :
                dist = get_distance(coo1,coo2)
                if dist <= contact[2]:
                    residue_list.append(True)
                else:
                    residue_list.append(False)
            else:
                print "Error in progressing"
    hits =0
    for entry in residue_list:
        if entry:
            hits +=1

    return contacts_in_smotif, float(hits)/float(contacts_in_smotif)*100.00

def S2ContactPredicition(s1_def, s2_def, smotif, exp_data, threshold):

    """

    :param s1_def:
    :param s2_def:
    :param smotif:
    :param threshold:
    :return:
    """
    #     print s1_def, s2_def ['strand', 9, 4, 4, 89, 97] ['strand', 9, 4, 4, 102, 110]

    ss1_list = range(s1_def[4], s1_def[5]+1)
    ss2_list = range(s2_def[4], s2_def[5]+1)

    smotif_ss1 = range(int(smotif[0][0][1]), int(smotif[0][0][2])+1 )
    smotif_ss2 = range(int(smotif[0][0][3]), int(smotif[0][0][4])+1 )

    # print ss1_list, ss2_list
    # print smotif_ss1, smotif_ss2
    # print smotif[0][0]

    contacts = exp_data['contacts']

    residue_list = []
    contacts_in_smotif = 0
    for contact in contacts:
        if contact[0] in ss1_list and contact[1] in ss2_list:
            contacts_in_smotif +=1
            res1_index = ss1_list.index(contact[0])
            res2_index = ss2_list.index(contact[1])
            #print contact, res1_index, res2_index

            smotif_res1 = smotif_ss1[res1_index]
            smotif_res2 = smotif_ss2[res2_index]
            #print smotif_res1, smotif_res2

            for res in smotif[0][1]:
                if res[2] == 'CA' and res[0] == smotif_res1 :
                   coo1 = [res[3], res[4], res[5]]

            for res in smotif[0][2]:
                if res[2] == 'CA' and res[0] == smotif_res2 :
                   coo2 = [res[3], res[4], res[5]]
            if coo1 and coo2 :
                dist = get_distance(coo1,coo2)
                if dist <= contact[2]:
                    residue_list.append(True)
                else:
                    residue_list.append(False)
            else:
                print "Error in progressing"
    hits =0
    for entry in residue_list:
        if entry:
            hits +=1

    return contacts_in_smotif, float(hits)/float(contacts_in_smotif)*100.00
