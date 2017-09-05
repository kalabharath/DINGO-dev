from filters.constraints.looplengthConstraint import get_dist

def getHCoorMatrix(ss_list, smotif):
    noe_matrix = {}
    count = 0
    for i in range(1, len(smotif), 5):
        coo1 = [smotif[i][3], smotif[i][4], smotif[i][5]]
        noe_matrix[ss_list[count]] = coo1
        count += 1
    return noe_matrix


def S1NOEprob(s1_def, s2_def, smotif, exp_data):
    """

    :param s1_def:
    :param s2_def:
    :param smotif:
    :param exp_data:
    :return:
    """

    noe_data = exp_data['noe_data']
    noes = noe_data[0]
    total_noes = noe_data[1]
    noes_found = 0.0
    smotif_noes = 0.0
    ss1_list = range(s1_def[4], s1_def[5] + 1)
    ss2_list = range(s2_def[4], s2_def[5] + 1)

    #smotif_ss1 = range(int(smotif[0][1]), int(smotif[0][2]) + 1)
    #smotif_ss2 = range(int(smotif[0][3]), int(smotif[0][4]) + 1)

    coor_matrix = getHCoorMatrix(ss1_list,  smotif[1])
    coor2_matrix = getHCoorMatrix(ss2_list, smotif[2])
    coor_matrix.update(coor2_matrix)
    resi = coor_matrix.keys()

    for noe in noes:
        res1, atm1, res2, atom2, noe_cutoff, tol = noe[0], noe[1], noe[2], noe[3], noe[4], noe[5]
        if res1 in resi and res2 in resi:
            smotif_noes += 1.0
            coo1 = coor_matrix[res1]
            coo2 = coor_matrix[res2]
            dist = get_dist(coo1, coo2)
            if dist <= noe_cutoff + tol:
                noes_found += 1.0

    if smotif_noes > 0:
        local_prob = (noes_found / smotif_noes)
    else:
        local_prob = 0.0
    #return (noes_found / total_noes), local_prob, noes_found
    return local_prob, local_prob, noes_found

def getSxCoorMatrix(coor_array, native_sse):

    resi ={}
    native_sse_range = range(native_sse[4], native_sse[5] + 1)
    count_array = 0
    for i in range(1, len(coor_array[0]), 5):
        x = (coor_array[0][i])
        y = (coor_array[1][i])
        z = (coor_array[2][i])
        resi[native_sse_range[count_array]] = [x, y, z]
        count_array += 1
    return resi

def SxNOEprob(transformed_coors, native_sse_order, current_ss, exp_data):
    """

    :param transformed_coors:
    :param native_sse_order:
    :param current_ss:
    :param exp_data:
    :return:
    """
    import copy
    sse_coors = copy.deepcopy(transformed_coors)

    noe_data = exp_data['noe_data']
    noes = noe_data[0]
    total_noes = noe_data[1]
    noes_found = 0.0
    smotif_noes = 0.0
    coor_matrix = {}

    for i in range(0, len(sse_coors)):
        tcoor = getSxCoorMatrix(sse_coors[i], native_sse_order[i])
        coor_matrix.update(tcoor)

    resi = coor_matrix.keys()

    for noe in noes:
        res1, res2, noe_cutoff, tol = noe[0], noe[2], noe[4], noe[5]

        if (res1 in resi) and (res2 in resi):  # this make sure that the resi are in the SSEs
            smotif_noes += 1.0
            coo1 = coor_matrix[res1]
            coo2 = coor_matrix[res2]
            dist = get_dist(coo1, coo2)
            if dist <= noe_cutoff+tol:
                noes_found += 1.0

    if smotif_noes > 0:
        local_prob = (noes_found / smotif_noes)
    else:
        local_prob = 0.0

    #return (noes_found / total_noes), local_prob, noes_found
    return local_prob, local_prob, noes_found
