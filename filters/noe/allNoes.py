from filters.contacts.contacts_filter import get_distance


def s1NOEfit(s1_def, s2_def, smotif, exp_data):
    noe_cutoff = False
    noe_matrix = exp_data['noe_data']
    ss1_list = range(s1_def[4], s1_def[5] + 1)
    ss2_list = range(s2_def[4], s2_def[5] + 1)

    smotif_ss1 = range(int(smotif[0][1]), int(smotif[0][2]) + 1)
    smotif_ss2 = range(int(smotif[0][3]), int(smotif[0][4]) + 1)

    noes_found = []
    noes_total = []

    # Scoring intra-SSE NOEs
    coo1 = [999, 999, 999]
    coo2 = [0.0, 0.0, 0.0]

    for i in range(0, len(smotif_ss1) - 1):
        sres1 = smotif_ss1[i]
        nres1 = ss1_list[i]
        for j in range(i + 1, len(smotif_ss1)):
            sres2 = smotif_ss1[j]
            nres2 = ss1_list[j]
            if noe_matrix[nres1, nres2]:
                noe_cutoff = noe_matrix[nres1, nres2]
                for entry1 in smotif[1]:
                    if entry1[2] == 'H' and entry1[0] == sres1:
                        coo1 = [entry1[3], entry1[4], entry1[5]]
                    if entry1[2] == 'H' and entry1[0] == sres2:
                        coo2 = [entry1[3], entry1[4], entry1[5]]

                dist = get_distance(coo1, coo2)
                # print sres1, sres2, noe_cutoff, coo1, coo2, dist
                if noe_cutoff > 10000:
                    real_noe = noe_cutoff - 10000
                    # backmapping side chain noes to amides
                    if (real_noe - 4.0 <= dist <= real_noe + 4.0):
                        noes_found.append((sres1, sres2))
                        noes_total.append((sres1, sres2))
                    else:
                        noes_total.append((sres1, sres2))
                elif dist <= noe_cutoff:
                    noes_found.append((sres1, sres2))
                    noes_total.append((sres1, sres2))
                else:
                    noes_total.append((sres1, sres2))
    coo1 = [999, 999, 999]
    coo2 = [0.0, 0.0, 0.0]

    for i in range(0, len(smotif_ss2) - 1):
        sres1 = smotif_ss2[i]
        nres1 = ss2_list[i]
        for j in range(i + 1, len(smotif_ss2)):
            sres2 = smotif_ss2[j]
            nres2 = ss2_list[j]
            if noe_matrix[nres1, nres2]:
                noe_cutoff = noe_matrix[nres1, nres2]
                for entry1 in smotif[2]:
                    if entry1[2] == 'H' and entry1[0] == sres1:
                        coo1 = [entry1[3], entry1[4], entry1[5]]
                    if entry1[2] == 'H' and entry1[0] == sres2:
                        coo2 = [entry1[3], entry1[4], entry1[5]]

                dist = get_distance(coo1, coo2)
                # print sres1, sres2, noe_cutoff, coo1, coo2, dist
                if noe_cutoff > 10000:
                    real_noe = noe_cutoff - 10000
                    # backmapping side chain noes to amides
                    if (real_noe - 4.0 <= dist <= real_noe + 4.0):
                        noes_found.append((sres1, sres2))
                        noes_total.append((sres1, sres2))
                    else:
                        noes_total.append((sres1, sres2))
                elif dist <= noe_cutoff:
                    noes_found.append((sres1, sres2))
                    noes_total.append((sres1, sres2))
                else:
                    noes_total.append((sres1, sres2))
    # end of Scoring intra-SSE NOEs


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
                                else:
                                    noes_total.append((res1, res2))
                            elif dist <= noe_cutoff:
                                noes_found.append((res1, res2))
                                noes_total.append((res1, res2))
                            else:
                                noes_total.append((res1, res2))
    if len(noes_found) == 0:
        return 0.00

    # fmeasure = calcFmeasure(noes_found, noes_total)
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
    noe_cutoff = False
    sse_satisfied = False
    sse_coors = copy.deepcopy(transformed_coors)
    noe_matrix = exp_data['noe_data']
    noes_found = []
    noes_total = []
    for i in range(0, len(sse_coors) - 1):
        # wtf am i doing here !
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
                    else:
                        noe_cutoff = False

                    if noe_cutoff:
                        ca_res2 = [ca2[0][ss2_list.index(res2)], ca2[1][ss2_list.index(res2)],
                                   ca2[2][ss2_list.index(res2)]]
                        dist = get_dist(ca_res1, ca_res2)

                        if noe_cutoff > 10000:
                            real_noe = noe_cutoff - 10000
                            # backmapping side chain noes to amides

                            if (real_noe - 4.0 <= dist <= real_noe + 4.0):
                                noes_found.append((res1, res2))
                                noes_total.append((res1, res2))
                            else:
                                noes_total.append((res1, res2))

                        elif dist <= noe_cutoff:
                            sse_satisfied = True
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
