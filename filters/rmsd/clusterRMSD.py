
from qcp import *


def combine_arrays(array1, array2):

    tarray = copy.deepcopy(array1)
    for i in range(0, len(tarray)):
        tarray[i] = tarray[i] + array2[i]
    return tarray


def copy_equivalent(frag, eq):

    #print "Wtf is worng with these arrays", len(frag[0]), len(eq)
    factor = len(eq)
    x, y, z = [None] * factor, [None] * factor, [None] * factor
    index = 0

    for i in eq:

        x[index] = (frag[0][i])
        y[index] = (frag[1][i])
        z[index] = (frag[2][i])
        index += 1
    return [x, y, z]


def parseCAequivalent(coor1, coor2, eq1, eq2):

    frag_aca = getCAcoo(coor1[0])

    for i in range(1, len(coor1)):
        tfraga = getCAcoo(coor1[i])
        frag_aca = combine_arrays(frag_aca, tfraga)

    frag_bca = getCAcoo(coor2[0])

    for i in range(1, len(coor2)):
        tfragb = getCAcoo(coor2[i])
        frag_bca = combine_arrays(frag_bca, tfragb)

    frag_aca = copy_equivalent(frag_aca, eq1)
    frag_bca = copy_equivalent(frag_bca, eq2)

    return frag_aca, frag_bca


def clusterRMSD(coor1, coor2, eq1, eq2):

    try:
        frag_aca, frag_bca = parseCAequivalent(coor1, coor2, eq1, eq2)
    except IndexError:
        return 999.999

    fraglen = len(frag_aca[0])
    xyz1 = qcprot.MakeDMatrix(3, fraglen)
    xyz2 = qcprot.MakeDMatrix(3, fraglen)

    frag_aca, a_cen = centerCoo(frag_aca)
    frag_bca, b_cen = centerCoo(frag_bca)

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

    qcprot.FreeDMatrix(xyz1)
    qcprot.FreeDMatrix(xyz2)
    qcprot.FreeDArray(rot)

    return rmsd