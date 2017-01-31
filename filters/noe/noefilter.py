import copy

from filters.constraints.looplengthConstraint import get_dist
from filters.contacts.contacts_filter import get_distance
from filters.contacts.evfoldContacts import calcFmeasure


def calcPrecision(gbar, ganoe):
    tp, fp, fn = 0.0, 0.0, 0.0

    for entry in ganoe:
        if entry in gbar:
            tp += 1
        else:
            fn += 1
    for entry in gbar:
        if entry in ganoe:
            pass
        else:
            fp += 1
    if tp:
        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        return precision
        # return (2 * precision * recall) / (precision + recall)
    else:
        return False


def s1NOEfit(s1_def, s2_def, smotif, exp_data):
    noe_cutoff = False
    noe_matrix = exp_data['noe_data']
    ss1_list = range(s1_def[4], s1_def[5] + 1)
    ss2_list = range(s2_def[4], s2_def[5] + 1)

    smotif_ss1 = range(int(smotif[0][1]), int(smotif[0][2]) + 1)
    smotif_ss2 = range(int(smotif[0][3]), int(smotif[0][4]) + 1)

    noes_found = []
    noes_total = []
    for res in smotif_ss1:
        for entry1 in smotif[1]:
            if entry1[2] == 'H' and entry1[0] == res:
                coo1 = [entry1[3], entry1[4], entry1[5]]
                for entry2 in smotif[2]:
                    if entry2[2] == 'H':
                        coo2 = [entry2[3], entry2[4], entry2[5]]
                        res1 = ss1_list[smotif_ss1.index(entry1[0])]
                        res2 = ss2_list[smotif_ss2.index(entry2[0])]
                        if noe_matrix[res1, res2]:
                            noe_cutoff = noe_matrix[res1, res2]
                        elif noe_matrix[res2, res1]:
                            noe_cutoff = noe_matrix[res2, res1]
                        else:
                            noe_cutoff = False

                        if noe_cutoff:
                            dist = get_distance(coo1, coo2)
                            if noe_cutoff > 10000:
                                real_noe = noe_cutoff - 10000
                                # backmapping side chain noes to amides
                                if (real_noe - 4.0 <= dist <= real_noe + 4.0):

                                    noes_found.append((res1, res2))
                                    noes_total.append((res1, res2))
                            elif dist <= noe_cutoff:
                                noes_found.append((res1, res2))
                                noes_total.append((res1, res2))
                            else:
                                noes_total.append((res1, res2))
    if len(noes_found) == 0:
            return 0.00


    """
    noes_total = []
    for entry1 in ss1_list:
        for entry2 in ss2_list:
            if noe_matrix[entry1, entry2]:
                noes_total.append((entry1, entry2))

    fmeasure = calcFmeasure(noes_found, noes_total)
    # precision = calcPrecision(noes_found, total_noes)
    # print len(total_noes), total_noes
    # print len(noes_found), noes_found
    """
    #fmeasure = calcFmeasure(noes_found, noes_total)
    fmeasure = (float(len(noes_found)) / float(len(noes_total)))
    return fmeasure


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


def s2NOEfit(transformed_coors, native_sse_order, exp_data):
    import warnings
    """
    Depr
    :param transformed_coors:
    :param native_sse_order:
    :param exp_data:
    :return:
    """
    warnings.warn("deprecated", DeprecationWarning)

    noe_cutoff = False

    sse_coors = copy.deepcopy(transformed_coors)
    noe_matrix = exp_data['noe_data']
    noes_found = []
    noes_total = []
    for i in range(0, len(sse_coors) - 1):
        res_c, ca1 = getNHandresi(sse_coors[i])
        res_n, ca2 = getNHandresi(sse_coors[i + 1])
        ss1_list = range(native_sse_order[i][4], native_sse_order[i][5] + 1)
        ss2_list = range(native_sse_order[i + 1][4], native_sse_order[i + 1][5] + 1)
        try:
            for res1 in ss1_list:
                # print ss1_list, res1, res_c, res_c[ss1_list.index(res1)]
                # print ss1_list.index(res1)
                ca_res1 = [ca1[0][ss1_list.index(res1)], ca1[1][ss1_list.index(res1)], ca1[2][ss1_list.index(res1)]]
                for res2 in ss2_list:
                    if noe_matrix[res1, res2]:
                        noe_cutoff = noe_matrix[res1, res2]
                    elif noe_matrix[res2, res1]:
                        noe_cutoff = noe_matrix[res2, res1]
                    else:
                        noe_cutoff = False
                    if noe_cutoff:
                        ca_res2 = [ca2[0][ss2_list.index(res2)], ca2[1][ss2_list.index(res2)],
                                   ca2[2][ss2_list.index(res2)]]
                        dist = get_dist(ca_res1, ca_res2)
                        if dist < noe_cutoff:
                            noes_found.append((res1, res2))
                            noes_total.append((res1, res2))
                        else:
                            noes_total.append((res1, res2))
        except:
            return []
    if len(noes_found) == 0:
        return []

    fmeasure = calcFmeasure(noes_found, noes_total)

    return fmeasure


def s3NOEfit(transformed_coors, native_sse_order, current_ss, exp_data):
    noe_cutoff = False

    sse_coors = copy.deepcopy(transformed_coors)
    noe_matrix = exp_data['noe_data']
    noes_found = []
    noes_total = []


    # current_ss ['strand', 13, 4, 10, 36, 48]

    cresi = range(current_ss[4], current_ss[5] + 1)

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
        if ss1_list == cresi or ss2_list == cresi:
            sse_satisfied = False
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

                        if noe_cutoff >10000:
                            real_noe = noe_cutoff-10000
                            #backmapping side chain noes to amides

                            if (real_noe - 5.0 <= dist <= real_noe + 5.0):

                                noes_found.append((res1, res2))
                                noes_total.append((res1, res2))

                        elif dist <= noe_cutoff :
                            sse_satisfied = True
                            noes_found.append((res1, res2))
                            noes_total.append((res1, res2))
                        else:
                            noes_total.append((res1, res2))

    #if len(noes_found) == 0 or not sse_satisfied:
    if len(noes_found) == 0:
        return 0.00
    #fmeasure = calcFmeasure(noes_found, noes_total)
    fmeasure = (float(len(noes_found)) / float(len(noes_total)))

    return fmeasure
