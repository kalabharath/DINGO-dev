"""
[['strand', 8, 6, 7, 14], ['strand', 7, 2, 17, 23], ['strand', 9, 4, 28, 36], ['strand', 7, 5, 42, 48],
 ['strand', 5, 2, 51, 55], ['helix', 6, 6, 62, 67], ['strand', 5, 5, 75, 79], ['strand', 7, 4, 84, 90]]
{0: [['strand', 8, 6, 2, 7, 14], ['strand', 6, 8, 2, 9, 14], ['strand', 7, 7, 2, 8, 14], ['strand', 10, 5, 1, 6, 15], ['strand', 9, 5, 2, 6, 14], ['strand', 10, 4, 2, 5, 14], ['strand', 6, 6, 4, 7, 12], ['strand', 7, 6, 3, 7, 13], ['strand', 6, 7, 3, 8, 13], ['strand', 9, 6, 1, 7, 15], ['strand', 10, 6, 0, 7, 16]], 1: [['strand', 7, 2, 4, 17, 23], ['strand', 5, 4, 4, 19, 23], ['strand', 6, 3, 4, 18, 23], ['strand', 9, 1, 3, 16, 24], ['strand', 8, 1, 4, 16, 23], ['strand', 9, 0, 4, 15, 23], ['strand', 5, 2, 6, 17, 21], ['strand', 6, 2, 5, 17, 22], ['strand', 5, 3, 5, 18, 22], ['strand', 8, 2, 3, 17, 24], ['strand', 9, 2, 2, 17, 25]], 2: [['strand', 9, 4, 5, 28, 36], ['strand', 7, 6, 5, 30, 36], ['strand', 8, 5, 5, 29, 36], ['strand', 11, 3, 4, 27, 37], ['strand', 10, 3, 5, 27, 36], ['strand', 11, 2, 5, 26, 36], ['strand', 7, 4, 7, 28, 34], ['strand', 8, 4, 6, 28, 35], ['strand', 7, 5, 6, 29, 35], ['strand', 10, 4, 4, 28, 37], ['strand', 11, 4, 3, 28, 38]], 3: [['strand', 7, 5, 2, 42, 48], ['strand', 5, 7, 2, 44, 48], ['strand', 6, 6, 2, 43, 48], ['strand', 9, 4, 1, 41, 49], ['strand', 8, 4, 2, 41, 48], ['strand', 9, 3, 2, 40, 48], ['strand', 5, 5, 4, 42, 46], ['strand', 6, 5, 3, 42, 47], ['strand', 5, 6, 3, 43, 47], ['strand', 8, 5, 1, 42, 49], ['strand', 9, 5, 0, 42, 50]], 4: [['strand', 5, 2, 6, 51, 55], ['strand', 7, 1, 5, 50, 56], ['strand', 6, 1, 6, 50, 55], ['strand', 7, 0, 6, 49, 55], ['strand', 6, 2, 5, 51, 56], ['strand', 7, 2, 4, 51, 57]], 5: [['helix', 6, 6, 5, 62, 67], ['helix', 5, 7, 5, 63, 67], ['helix', 8, 5, 4, 61, 68], ['helix', 7, 5, 5, 61, 67], ['helix', 8, 4, 5, 60, 67], ['helix', 5, 6, 6, 62, 66], ['helix', 7, 6, 4, 62, 68], ['helix', 8, 6, 3, 62, 69]], 6: [['strand', 5, 5, 4, 75, 79], ['strand', 7, 4, 3, 74, 80], ['strand', 6, 4, 4, 74, 79], ['strand', 7, 3, 4, 73, 79], ['strand', 6, 5, 3, 75, 80], ['strand', 7, 5, 2, 75, 81]], 7: [['strand', 7, 4, 24, 84, 90], ['strand', 5, 6, 24, 86, 90], ['strand', 6, 5, 24, 85, 90], ['strand', 9, 3, 23, 83, 91], ['strand', 8, 3, 24, 83, 90], ['strand', 9, 2, 24, 82, 90], ['strand', 5, 4, 26, 84, 88], ['strand', 6, 4, 25, 84, 89], ['strand', 5, 5, 25, 85, 89], ['strand', 8, 4, 23, 84, 91], ['strand', 9, 4, 22, 84, 92]]}


"""


def genSSDef(ss_seq):
    """
	returns an array of smotifs derived from the ss_seq with extended
	functionality ie information about left_loop length
	return [[ss_type,len_ss,l_loop,start,end]['strand', 5, 4, 5, 9]]
	:param ss_seq:
	:return:
	"""
    helix, strand, loop, l_loop = 1, 1, 1, 0
    return_array = []
    seq = []
    for i in range(0, len(ss_seq)):
        if i == 0:
            continue
        else:
            if (ss_seq[i] == 'C') or (ss_seq[i] == 'L'):
                l_loop += 1
            if ss_seq[i] == ss_seq[i - 1]:
                if ss_seq[i] == 'H':
                    seq.append(i)
                    helix += 1
                if ss_seq[i] == 'E':
                    seq.append(i)
                    strand += 1
                if (ss_seq[i] == 'C') or (ss_seq[i] == 'L'):
                    loop += 1
            else:
                if (helix > 4) and ((strand == 1)):
                    end = i
                    return_array.append(['helix', helix, l_loop, seq[0], end])
                    seq = []
                    l_loop = 0
                elif (strand > 4) and (helix == 1):
                    end = i
                    return_array.append(['strand', strand, l_loop, seq[0], end])
                    seq = []
                    l_loop = 0
                elif (loop > 4) and ((strand == 1) and (helix == 1)):
                    end = i
                    seq = []
                helix, strand, loop = 1, 1, 1
                seq = []
    return return_array, l_loop + 1


def getCompleteSSDef(ss_seq):
    exn_ss, last_loop = genSSDef(ss_seq)
    import copy
    com_ss = []
    for i in range(0, len(exn_ss) - 1):

        current_def = exn_ss[i]
        next_def = exn_ss[i + 1]
        txx = copy.copy(current_def)
        txx.insert(3, next_def[2])
        com_ss.append(txx)
        if i == len(exn_ss) - 2:
            txx = copy.copy(next_def)
            txx.insert(3, last_loop)
            com_ss.append(txx)
    return com_ss, exn_ss


def genPermutations(exn_ss):
    ss_combi = {}
    strict_2 = [(0, 0), (-2, 0), (-2, 1), (-2, 2), (-1, -1), (-1, 0), (-1, 1), (-1, 2), (0, -2), (0, -1), (0, 1),
                (0, 2), (1, -2), (1, -1), (1, 0), (1, 1), (2, -2), (2, -1), (2, 0)]

    for  k in range(0, len(exn_ss)):
        ss_type = exn_ss[k][0]
        len_ss = exn_ss[k][1]
        l_loop = exn_ss[k][2]
        r_loop = exn_ss[k][3]
        start = exn_ss[k][4]
        end = exn_ss[k][5]
        #print exn_ss[k]
        for npr in strict_2:
            nter = npr[0]
            cter = npr[1]
            # print nter, cter
            tl_loop = l_loop
            tr_loop = r_loop
            tstart = start
            tend = end
            tlen_ss = len_ss

            if nter > 0:
                tl_loop = l_loop - nter
                if tl_loop >= 1:
                    tstart = tstart - nter
                else:
                    continue
            elif nter < 0:
                tl_loop = l_loop + ((-1) * nter)
                if tl_loop >= 1:
                    tstart = tstart + ((-1) * nter)
                else:
                    continue
            else:
                pass

            if cter > 0:
                tr_loop = r_loop - cter
                if tr_loop >= 1:
                    tend = tend + cter
                else:
                    continue
            elif cter < 0:
                tr_loop = r_loop + ((-1) * cter)
                if tr_loop >= 1:
                    tend = tend - ((-1) * cter)
                else:
                    continue
            else:
                pass

            tlen_ss = len_ss + nter + cter
            arr_0 = [ss_type, tlen_ss, tl_loop, tr_loop, tstart, tend]

            ss_combi.setdefault(k, []).append(arr_0)


    return ss_combi


def genSSPermutations(ss_seq):
    exn_ss, ss_def = getCompleteSSDef(ss_seq)
    ss_combi = genPermutations(exn_ss)

    return ss_def, ss_combi
