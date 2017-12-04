
import copy
import qcprot
from qcp import *


def combine_arrays(array1, array2):
    tarray = copy.deepcopy(array1)
    for i in range(0, len(tarray)):
        tarray[i] = tarray[i] + array2[i]
    return tarray


def refineRMSD(smotif_coors, pair, csmotif, rmsd_cutoff):

    frag_a1 = copy.deepcopy(smotif_coors[pair[0]])
    frag_a2 = copy.deepcopy(smotif_coors[pair[1]])
    frag_a = combine_arrays(frag_a1, frag_a2)

    frag_b1 = getcoo(csmotif[1])
    frag_b2 = getcoo(csmotif[2])
    frag_b = combine_arrays(frag_b1, frag_b2)

    frag_a, a_cen = centerCoo(frag_a)
    frag_b, b_cen = centerCoo(frag_b)

    frag_aca = getCAcoo(frag_a)
    frag_bca = getCAcoo(frag_b)

    fraglen = len(frag_aca[0])

    xyz1 = qcprot.MakeDMatrix(3, fraglen)
    xyz2 = qcprot.MakeDMatrix(3, fraglen)

    for i in range(0, fraglen):
        qcprot.SetDArray(0, i, xyz1, frag_aca[0][i])
        qcprot.SetDArray(1, i, xyz1, frag_aca[1][i])
        qcprot.SetDArray(2, i, xyz1, frag_aca[2][i])
    for i in range(0, fraglen):
        qcprot.SetDArray(0, i, xyz2, frag_bca[0][i])
        qcprot.SetDArray(1, i, xyz2, frag_bca[1][i])
        qcprot.SetDArray(2, i, xyz2, frag_bca[2][i])

    rot = qcprot.MakeDvector(9)

    # *********
    rmsd = qcprot.CalcRMSDRotationalMatrix(xyz1, xyz2, fraglen, rot)
    # *********

    if rmsd > rmsd_cutoff:
        # exit without further computation
        # free memory
        qcprot.FreeDMatrix(xyz1)
        qcprot.FreeDMatrix(xyz2)
        qcprot.FreeDArray(rot)
        return [], rmsd


    rotmat = [None] * 9
    for i in range(0, 9):
        rotmat[i] = qcprot.GetDvector(i, rot)

    cm_frag_b1 = translateCM(frag_b1, b_cen)
    cm_frag_b2 = translateCM(frag_b2, b_cen)

    rot_frag_b1 = applyRot(cm_frag_b1, rotmat)
    rot_frag_b2 = applyRot(cm_frag_b2, rotmat)

    trans_frag_b1 = applyTranslation(rot_frag_b1, a_cen)
    trans_frag_b2 = applyTranslation(rot_frag_b2, a_cen)

    smotif_coors.pop(pair[0])
    smotif_coors.insert(pair[0], trans_frag_b1)

    smotif_coors.pop(pair[1])
    smotif_coors.insert(pair[1], trans_frag_b2)


    qcprot.FreeDMatrix(xyz1)
    qcprot.FreeDMatrix(xyz2)
    qcprot.FreeDArray(rot)


    return smotif_coors, rmsd


def kClahsesRefined(coo_arrays, sse_ordered, sse_pairs):

    sse1 = sse_ordered[sse_pairs[0]]
    sse2 = sse_ordered[sse_pairs[1]]

    if kClashes(coo_arrays, sse_ordered, sse1):
        if kClashes(coo_arrays, sse_ordered, sse2):
            return True
        else:
            return False

    else:
        return False




def kClashesOLD(coo_arrays, sse_ordered):
    """
    Known minimum distances between various atoms between the SSEs of the smotif
    the first entry is mean and second entry is SD, computed on a sample of 100 lowest distances observed over single
    digit smotifs that have the largest number of entries.
    :param coo_arrays:
    :param cdist:
    :param sse_ordered:
    :return:
    """

    # [['strand', 15, 10, 3, 59, 73], ['strand', 15, 3, 13, 77, 91], ['strand', 14, 11, 4, 103, 116]]
    # atom_array = ['N', 'H', 'CA', 'C', 'O']

    for i in range(0, len(coo_arrays) - 1):
        # Compute clashes for amide hydrogen pairs
        atom_type = 1
        sse1 = getXcoo(coo_arrays[i], atom_type)
        for p in range(i + 1, len(coo_arrays)):
            kdist = getKdist([sse_ordered[i], sse_ordered[p]], atom_type)
            sse2 = getXcoo(coo_arrays[p], atom_type)
            for j in range(0, len(sse1[0])):
                for k in range(0, len(sse2[0])):
                    dist = get_dist([sse1[0][j], sse1[1][j], sse1[2][j]], [sse2[0][k], sse2[1][k], sse2[2][k]])
                    dist = round(dist, 2)
                    if dist < kdist:
                        return False

    for i in range(0, len(coo_arrays) - 1):
        # Compute clashes for carbonyl carbon pairs
        atom_type = 3
        sse1 = getXcoo(coo_arrays[i], atom_type)
        for p in range(i + 1, len(coo_arrays)):
            kdist = getKdist([sse_ordered[i], sse_ordered[p]], atom_type)
            sse2 = getXcoo(coo_arrays[p], atom_type)
            for j in range(0, len(sse1[0])):
                for k in range(0, len(sse2[0])):
                    dist = get_dist([sse1[0][j], sse1[1][j], sse1[2][j]], [sse2[0][k], sse2[1][k], sse2[2][k]])
                    dist = round(dist, 2)
                    if dist < kdist:
                        return False
    return True