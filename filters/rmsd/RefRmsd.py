import copy

import qcp
import qcprot


def getRefCA(ref_ca, s1_def, s2_def):
    x, y, z = [], [], []

    for res in s1_def:
        x.append(ref_ca[int(res) - 1][0])
        y.append(ref_ca[int(res) - 1][1])
        z.append(ref_ca[int(res) - 1][2])

    for res in s2_def:
        x.append(ref_ca[int(res) - 1][0])
        y.append(ref_ca[int(res) - 1][1])
        z.append(ref_ca[int(res) - 1][2])

    return [x,y,z]

def concatArrays(frag_aca, frag_bca):
    xx = frag_aca[0] + frag_bca[0]
    yy = frag_aca[1] + frag_bca[1]
    zz = frag_aca[2] + frag_bca[2]
    return [xx, yy, zz]

def calcRefRMSD(ref_ca, s1_def, s2_def, smotif, rmsd_cutoff):
    # print 'ref_ca', ref_ca
    # print  'S1_def', s1_def
    # print  'S2_def', s2_def
    # print  'Smotif_data', smotif

    ss1_list = range(s1_def[4], s1_def[5] + 1)
    ss2_list = range(s2_def[4], s2_def[5] + 1)


    smotif_ss1 = range(int(smotif[0][1]), int(smotif[0][2]) + 1)
    smotif_ss2 = range(int(smotif[0][3]), int(smotif[0][4]) + 1)

    # print ss1_list, ss2_list
    #print smotif_ss1, smotif_ss2

    ref_coo = getRefCA(ref_ca, ss1_list, ss2_list)

    ref_coo, ref_coo_cen = qcp.centerCoo(ref_coo)



    frag_a = qcp.getcoo(smotif[1])
    frag_b = qcp.getcoo(smotif[2])

    frag_a, a_cen = qcp.centerCoo(frag_a)
    frag_b, b_cen = qcp.centerCoo(frag_b)

    frag_aca = qcp.getCAcoo(frag_a)
    frag_bca = qcp.getCAcoo(frag_b)
    smotif_coo = concatArrays(frag_aca, frag_bca)

    fraglen = len(smotif_coo[0])
    xyz1 = qcprot.MakeDMatrix(3, fraglen)
    xyz2 = qcprot.MakeDMatrix(3, fraglen)

    for i in range(0, fraglen):
        qcprot.SetDArray(0, i, xyz1, smotif_coo[0][i])
        qcprot.SetDArray(1, i, xyz1, smotif_coo[1][i])
        qcprot.SetDArray(2, i, xyz1, smotif_coo[2][i])
    for i in range(0, fraglen):
        qcprot.SetDArray(0, i, xyz2, ref_coo[0][i])
        qcprot.SetDArray(1, i, xyz2, ref_coo[1][i])
        qcprot.SetDArray(2, i, xyz2, ref_coo[2][i])

    rot = qcprot.MakeDvector(9)
    # *********
    rmsd = qcprot.CalcRMSDRotationalMatrix(xyz1, xyz2, fraglen, rot)
    # *********
    qcprot.FreeDMatrix(xyz1)
    qcprot.FreeDMatrix(xyz2)
    qcprot.FreeDArray(rot)

    if rmsd < rmsd_cutoff:
        return rmsd
    else:
        return False

def getRefCA2(ref_ca, sse_ordered):
    x, y, z = [], [], []
    for sse in sse_ordered:
        sse_def = range(sse[4], sse[5] + 1)
        for res in sse_def:
            x.append(ref_ca[int(res) - 1][0])
            y.append(ref_ca[int(res) - 1][1])
            z.append(ref_ca[int(res) - 1][2])

    return [x, y, z]


def calcRefRMSD2(ref_ca, sse_ordered, transformed_coos):

    ref_coo = getRefCA2(ref_ca, sse_ordered)
    ref_coo, ref_coo_cen = qcp.centerCoo(ref_coo)

    sse_coors = copy.copy(transformed_coos)

    smotif_coo = [[],[],[]]
    for sse in sse_coors:
        #frag, a_cen = qcp.centerCoo(sse)
        frag_ca = qcp.getCAcoo(sse)
        smotif_coo = concatArrays(smotif_coo, frag_ca)

    smotif_coo, xx_cen = qcp.centerCoo(smotif_coo)
    # print len(smotif_coo[0]), len(ref_coo[0])


    fraglen = len(smotif_coo[0])
    xyz1 = qcprot.MakeDMatrix(3, fraglen)
    xyz2 = qcprot.MakeDMatrix(3, fraglen)

    for i in range(0, fraglen):
        qcprot.SetDArray(0, i, xyz1, smotif_coo[0][i])
        qcprot.SetDArray(1, i, xyz1, smotif_coo[1][i])
        qcprot.SetDArray(2, i, xyz1, smotif_coo[2][i])
    for i in range(0, fraglen):
        qcprot.SetDArray(0, i, xyz2, ref_coo[0][i])
        qcprot.SetDArray(1, i, xyz2, ref_coo[1][i])
        qcprot.SetDArray(2, i, xyz2, ref_coo[2][i])

    rot = qcprot.MakeDvector(9)
    # *********
    rmsd = qcprot.CalcRMSDRotationalMatrix(xyz1, xyz2, fraglen, rot)
    # *********
    qcprot.FreeDMatrix(xyz1)
    qcprot.FreeDMatrix(xyz2)
    qcprot.FreeDArray(rot)

    return rmsd
