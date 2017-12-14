#!/usr/bin/env python

"""
Project_Name: constraints/looplengthConstraint, File_name: looplengthConstraint.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 27/11/15 , Time:2:28 PM
"""
from filters.rmsd.qcp import getCAcoo

def get_dist(r1, r2):
    import math
    x1, y1, z1 = r1[0], r1[1], r1[2]
    x2, y2, z2 = r2[0], r2[1], r2[2]
    return math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1))


def loopConstraint(coo_arrays, sseorder, direction, smotif_def):
    """
    Determine whether the assembled Smotif is too far away
    3.4 A * number of loop residues
    :param coo_arrays:
    :param sseorder:
    :param direction:
    :return:
    """

    nsh_dict = [0, 3.809, 3.137, 2.818, 2.482, 2.154, 1.928, 1.749, 1.67, 1.531, 1.428, 1.377, 1.282, 1.261, 1.203,
                1.135, 1.045, 1.004, 1.02, 0.977, 0.928, 0.865, 0.834, 0.811, 0.756, 0.761, 0.749, 0.777, 0.74, 0.655,
                0.648]
    nhs_dict = [0, 3.809, 3.137, 2.818, 2.482, 2.154, 1.928, 1.749, 1.67, 1.531, 1.428, 1.377, 1.282, 1.261, 1.203,
                1.135, 1.045, 1.004, 1.02, 0.977, 0.928, 0.865, 0.834, 0.811, 0.756, 0.761, 0.749, 0.777, 0.74, 0.655,
                0.648]
    nhh_dict = [0, 3.81, 3.036, 2.836, 2.511, 2.275, 2.178, 2.026, 1.876, 1.835, 1.669, 1.658, 1.666, 1.625, 1.53,
                1.445, 1.374, 1.292, 1.212, 1.164, 1.133, 1.049, 1.043, 1.074, 0.977, 0.965, 0.938, 0.868, 0.824, 0.805,
                0.788]
    nss_dict = [0, 3.81, 3.19, 1.846, 1.607, 1.274, 1.14, 1.139, 1.198, 1.177, 1.115, 1.029, 1.048, 0.935, 0.91, 0.908,
                0.85, 0.83, 0.852, 0.849, 0.761, 0.722, 0.742, 0.684, 0.677, 0.611, 0.587, 0.596, 0.565, 0.576, 0.532]
    hh_std = [0, 0.027, 0.284, 0.397, 0.441, 0.483, 0.499, 0.504, 0.537, 0.534, 0.538, 0.545, 0.507, 0.494, 0.468,
              0.447, 0.428, 0.439, 0.415, 0.432, 0.392, 0.382, 0.38, 0.401, 0.381, 0.38, 0.317, 0.328, 0.304, 0.318,
              0.273]
    ss_std = [0, 0.027, 0.313, 0.293, 0.469, 0.419, 0.474, 0.49, 0.505, 0.447, 0.501, 0.475, 0.479, 0.417, 0.451, 0.416,
              0.373, 0.395, 0.47, 0.418, 0.36, 0.349, 0.359, 0.312, 0.302, 0.281, 0.279, 0.264, 0.259, 0.346, 0.257]
    sh_std = [0, 0.067, 0.278, 0.361, 0.418, 0.45, 0.448, 0.455, 0.436, 0.452, 0.438, 0.416, 0.407, 0.402, 0.411, 0.405,
              0.381, 0.378, 0.373, 0.36, 0.372, 0.338, 0.322, 0.308, 0.285, 0.289, 0.296, 0.298, 0.294, 0.286, 0.208]
    hs_std = [0, 0.067, 0.278, 0.361, 0.418, 0.45, 0.448, 0.455, 0.436, 0.452, 0.438, 0.416, 0.407, 0.402, 0.411, 0.405,
              0.381, 0.378, 0.373, 0.36, 0.372, 0.338, 0.322, 0.308, 0.285, 0.289, 0.296, 0.298, 0.294, 0.286, 0.208]

    if direction == 'right':
        csse = sseorder[-1]
        psse = sseorder[-2]
        loop_length = csse[-2] - psse[-1]
        c_coo = getCAcoo(coo_arrays[-1])
        p_coo = getCAcoo(coo_arrays[-2])
        c_CA = [c_coo[0][0], c_coo[1][0], c_coo[2][0]]
        p_CA = [p_coo[0][-1], p_coo[1][-1], p_coo[2][-1]]

    else:
        csse = sseorder[0]
        psse = sseorder[1]
        loop_length = (psse[-2] - csse[-1])

        c_coo = getCAcoo(coo_arrays[0])
        p_coo = getCAcoo(coo_arrays[1])
        c_CA = [c_coo[0][-1], c_coo[1][-1], c_coo[2][-1]]
        p_CA = [p_coo[0][0], p_coo[1][0], p_coo[2][0]]

    dist = get_dist(c_CA, p_CA)

    if loop_length > 30.0 or loop_length == 0.0:
        return False

    Ndist = round(dist / float(loop_length), 2)

    stat_dist = 0
    stat_std = 0
    if smotif_def[0] == 'hh':
        stat_dist = nhh_dict[loop_length]
        stat_std = hh_std[loop_length]
    if smotif_def[0] == 'hs':
        stat_dist = nhs_dict[loop_length]
        stat_std = hs_std[loop_length]
    if smotif_def[0] == 'sh':
        stat_dist = nsh_dict[loop_length]
        stat_std = sh_std[loop_length]
    if smotif_def[0] == 'ss':
        stat_dist = nss_dict[loop_length]
        stat_std = ss_std[loop_length]

    if stat_dist - stat_std <= Ndist <= stat_dist + stat_std:
        return True
    else:
        return False


"""
    max_dist = loop_length * 3.4  # 3.8 Angstrom is the length between CA-CA in a linear polypeptide chain
    # however on a peptide triangular like arrangement the avg distance is 3.4 A

    if dist > max_dist:
        return False
    else:
        return True
"""
