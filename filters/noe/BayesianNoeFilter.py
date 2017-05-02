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

        res1, atm1, res2, atom2, noe_cutoff = noe[0], noe[1], noe[2], noe[3], noe[4]

        if res1 in resi and res2 in resi:
            smotif_noes += 1.0
            coo1 = coor_matrix[res1]
            coo2 = coor_matrix[res2]
            dist = get_dist(coo1, coo2)
            if dist <= noe_cutoff:
                noes_found += 1.0
    try:
        local_likelihood = (noes_found / smotif_noes)
    except:
        local_likelihood = 0.0

    return (noes_found / total_noes), local_likelihood, noes_found

def getNHandresi(frag):
    x, y, z = [], [], []
    resi = []
    for i in range(0, len(frag[0])):
        if frag[3][i] == 'H':
            x.append(frag[0][i])
            y.append(frag[1][i])
            z.append(frag[2][i])
            resi.append(frag[4][i])
    return resi, [x, y, z]

def sXGlobalNOEfit(transformed_coors, native_sse_order, current_ss, exp_data):
    import copy
    sse_coors = copy.deepcopy(transformed_coors)

    noe_data = exp_data['noe_data']
    noe_matrix = noe_data[0]
    total_noes = noe_data[1]

    noes_found = []
    noes_total = []

    # current_ss ['strand', 13, 4, 10, 36, 48]
    cresi = range(current_ss[4], current_ss[5] + 1)

    # Score Intra SSE NOEs first
    score_intra_sse = False
    for i in range(0, len(sse_coors)):
        ss1_list = range(native_sse_order[i][4], native_sse_order[i][5] + 1)
        res_c, ca1 = getNHandresi(sse_coors[i])
        try:
            for j in range(0, len(cresi) - 1):
                res1 = cresi[j]
                coo1 = [ca1[0][ss1_list.index(res1)], ca1[1][ss1_list.index(res1)], ca1[2][ss1_list.index(res1)]]
                for k in range(j + 1, len(cresi)):
                    res2 = cresi[k]
                    coo2 = [ca1[0][ss1_list.index(res2)], ca1[1][ss1_list.index(res2)],
                            ca1[2][ss1_list.index(res2)]]
                    noe_cutoff = noe_matrix[res1, res2]
                    if noe_cutoff:
                        dist = get_dist(coo1, coo2)
                        if noe_cutoff > 10000:
                            real_noe = noe_cutoff - 10000
                            # backmapping side chain noes to amides
                            if real_noe - 4.0 <= dist <= real_noe + 4.0:
                                noes_found.append((res1, res2))
                                noes_total.append((res1, res2))
                            else:
                                noes_total.append((res1, res2))
                        elif dist <= noe_cutoff:
                            noes_found.append((res1, res2))
                            noes_total.append((res1, res2))
                        else:
                            noes_total.append((res1, res2))
        except:
            continue

    import itertools
    sse_pairs = list(itertools.combinations(range(len(sse_coors)), 2))
    # [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (2, 3),\
    # (2, 4), (2, 5), (2, 6), (2, 7), (3, 4), (3, 5), (3, 6), (3, 7), (4, 5), (4, 6), (4, 7), (5, 6), (5, 7), (6, 7)]

    for ij in sse_pairs:
        i = ij[0]
        j = ij[1]
        res_c, ca1 = getNHandresi(sse_coors[i])
        res_n, ca2 = getNHandresi(sse_coors[j])
        ss1_list = range(native_sse_order[i][4], native_sse_order[i][5] + 1)
        ss2_list = range(native_sse_order[j][4], native_sse_order[j][5] + 1)

        for res1 in ss1_list:
            try:
                ca_res1 = [ca1[0][ss1_list.index(res1)], ca1[1][ss1_list.index(res1)], ca1[2][ss1_list.index(res1)]]
            except:
                continue
            for res2 in ss2_list:
                if noe_matrix[res1, res2]:
                    noe_cutoff = noe_matrix[res1, res2]
                elif noe_matrix[res2, res1]:
                    noe_cutoff = noe_matrix[res2, res1]
                else:
                    noe_cutoff = False
                if noe_cutoff:
                    # mo res1, res2, noe_matrix[res1, res2], noe_matrix[res2, res1], noe_cutoff
                    try:
                        ca_res2 = [ca2[0][ss2_list.index(res2)], ca2[1][ss2_list.index(res2)],
                                   ca2[2][ss2_list.index(res2)]]
                    except:
                        continue
                    dist = get_dist(ca_res1, ca_res2)

                    if noe_cutoff > 10000:
                        real_noe = noe_cutoff - 10000  # back mapping side chain noes to amides

                        if real_noe - 4.0 <= dist <= real_noe + 4.0:
                            noes_found.append((res1, res2))
                            noes_total.append((res1, res2))
                        else:
                            noes_total.append((res1, res2))

                    elif dist <= noe_cutoff:
                        noes_found.append((res1, res2))
                        noes_total.append((res1, res2))
                    else:
                        noes_total.append((res1, res2))

    # if len(noes_found) == 0 or not sse_satisfied:
    if len(noes_found) == 0:
        return 0.00, 0

    probability = (float(len(noes_found))) / float(total_noes)
    return probability, len(noes_found)
