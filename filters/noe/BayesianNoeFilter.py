from filters.constraints.looplengthConstraint import get_dist


def getCoorMatrix(ss_list, smotif_ss, smotif):
    noe_matrix = {}
    for i in range(0, len(smotif_ss)):
        actual_res = ss_list[i]
        smotif_res = smotif_ss[i]
        for entry in smotif:
            if entry[2] == 'H' and entry[0] == smotif_res:
                coo1 = [entry[3], entry[4], entry[5]]
                noe_matrix[actual_res] = coo1
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

    smotif_ss1 = range(int(smotif[0][1]), int(smotif[0][2]) + 1)
    smotif_ss2 = range(int(smotif[0][3]), int(smotif[0][4]) + 1)

    coor_matrix = getCoorMatrix(ss1_list, smotif_ss1, smotif[1])
    coor2_matrix = getCoorMatrix(ss2_list, smotif_ss2, smotif[2])
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
    return (noes_found / total_noes), local_prob, noes_found

def getSxCoorMatrix(coor_array, native_sse):

    resi ={}
    native_sse_range = range(native_sse[4], native_sse[5] + 1)
    count_array = 0
    for i in range(0, len(coor_array[0])):
        if coor_array[3][i] == 'H':
            x =  (coor_array[0][i])
            y =  (coor_array[1][i])
            z =  (coor_array[2][i])
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

        res1, atm1, res2, atom2, noe_cutoff, tol = noe[0], noe[1], noe[2], noe[3], noe[4], noe[5]

        if res1 in resi and res2 in resi:
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


    return (noes_found / total_noes), local_prob, noes_found
