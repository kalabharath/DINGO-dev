#!/usr/bin/env python

"""
Project_Name: constraints/looplengthConstraint, File_name: looplengthConstraint.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 27/11/15 , Time:2:28 PM
"""

def getCAcoo(frag):
    """
    :param frag:
    :return:
    """
    #print frag
    x, y, z = [], [], []
    for i in range(0, len(frag[0])):
        if frag[3][i] == 'CA':
            x.append(frag[0][i])
            y.append(frag[1][i])
            z.append(frag[2][i])
    return [x, y, z]

def get_dist(r1, r2):
    import math
    x1, y1, z1 = r1[0], r1[1], r1[2]
    x2, y2, z2 = r2[0], r2[1], r2[2]
    return math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))

def loopConstraint(coo_arrays, sseorder, direction):

    #  ['helix', 28, 3, 1, 78, 105] right
    # [['helix', 9, 4, 6, 5, 13], ['helix', 16, 4, 3, 20, 35], ['helix', 17, 3, 6, 39, 55], ['helix', 12, 5, 5, 63, 74], ['helix', 28, 4, 0, 79, 106]] right


    if direction == 'right':
        csse = sseorder[-1]
        psse = sseorder[-2]
        loop_length = csse[-2] - psse[-1]
        c_coo = getCAcoo(coo_arrays[-1])
        p_coo = getCAcoo(coo_arrays[-2])
        c_CA = [c_coo[0][0],c_coo[1][0],c_coo[2][0]]
        p_CA = [p_coo[0][-1],p_coo[1][-1],p_coo[2][-1]]


    else:
        csse = sseorder[0]
        psse = sseorder[1]
        loop_length = psse[-2] - csse[-1]
        c_coo = getCAcoo(coo_arrays[0])
        p_coo = getCAcoo(coo_arrays[1])
        c_CA = [c_coo[0][-1],c_coo[1][-1],c_coo[2][-1]]
        p_CA = [p_coo[0][0],p_coo[1][0],p_coo[2][0]]


    #print c_CA, p_CA

    dist = get_dist(c_CA, p_CA)

    max_dist = loop_length * 3.4 #3.8 Angstrom is the length between CA-CA in a linear polypeptide chain
    # however on a peptide triangular like arrangement the avg distance is 3.4 A

    if dist > max_dist:
        #print "Current distance: ", dist, max_dist
        return False
    else:
        return True

