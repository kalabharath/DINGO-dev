from filters.rmsd.qcp import *
from utility.io_util import readPickle

def processRBBC(cluster1):
    x, y, z = [None] * 5, [None] * 5, [None] * 5
    for i in range(0, 5):
        x[i] = cluster1[i][1]
        y[i] = cluster1[i][2]
        z[i] = cluster1[i][3]
    return [x, y, z]

def formatClusterCoo(cluster):
    tlen = len(cluster) *len(cluster[0])
    x, y, z, spin = [None] * tlen, [None] * tlen, [None] * tlen, [None] * tlen
    count = 0
    for i in range(0, len(cluster)):
        for j in range(0, len(cluster[i])):
            spin[count] = cluster[i][j][0]
            x[count] = cluster[i][j][1]
            y[count] = cluster[i][j][2]
            z[count] = cluster[i][j][3]
            count += 1

    return [x, y, z, spin]

def extractSpinCoo(cluster, spin, res_type):
    spin_indices = []
    repeat = []

    x, y, z, spin_type = [], [], [], []

    if res_type == 'I':
        ile = ['N', 'H', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1', 'HA', 'HB', '1HG1', '2HG1', '1HG2', '2HG2', '3HG2',
               '1HD1', '2HD1', '3HD1']
        repeat = len(ile)
        if spin =='HG2':
            spin_indices = [13, 14, 15]
        if spin == 'HD1':
            spin_indices = [16, 17, 18]
    if res_type == 'L':
        leu = ['N', 'H', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'HA', '1HB', '2HB', 'HG', '1HD1', '2HD1', '3HD1',
               '1HD2', '2HD2', '3HD2']
        repeat = len(leu)
        if spin == 'HD1':
            spin_indices = [13, 14, 15]
        if spin == 'HD2':
            spin_indices = [16, 17, 18]
        if spin == 'HD':
            spin_indices = [13, 14, 15, 16, 17, 18]
    if res_type == 'V':
        val = ['N', 'H', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'HA', 'HB', '1HG1', '2HG1', '3HG1', '1HG2', '2HG2', '3HG2']
        repeat = len(val)
        if spin == 'HG1':
            spin_indices = [10, 11, 12]
        if spin == 'HG2':
            spin_indices = [13, 14, 15]
        if spin == 'HG':
            spin_indices = [10, 11, 12, 13, 14, 15]

    if res_type == 'A':
        ala = ['N', 'H', 'CA', 'C', 'O', 'CB', 'HA', '1HB', '2HB', '3HB']
        repeat = len(ala)
        if spin == 'HB':
            spin_indices = [7, 8, 9]

    if spin_indices and repeat:
        for entry in spin_indices:
            for i in range(entry, len(cluster[0]), repeat):
                x.append(cluster[0][i])
                y.append(cluster[1][i])
                z.append(cluster[2][i])
                spin_type.append(cluster[3][i])

    return [x, y, z, spin_type]

def extend_array(extended, temp):
    if len(extended) == 0:
        return temp
    for i in range(0,len(temp)):
        extended[i] = extended[i]+ temp[i]

    return extended

def bbrmsd(bbc, rotamer_cluster, rmsd_cutoff, spin, res_type):


    fraga, a_cen = centerCoo(bbc)
    fraglen = 5
    all_spin_coors = []
    for cluster in rotamer_cluster:

        data = readPickle(cluster)
        rbbc = processRBBC(data[0])
        fragb, b_cen = centerCoo(rbbc)

        xyz1 = qcprot.MakeDMatrix(3, fraglen)
        xyz2 = qcprot.MakeDMatrix(3, fraglen)

        for i in range(0, fraglen):
            qcprot.SetDArray(0, i, xyz1, fraga[0][i])
            qcprot.SetDArray(1, i, xyz1, fraga[1][i])
            qcprot.SetDArray(2, i, xyz1, fraga[2][i])
        for i in range(0, fraglen):
            qcprot.SetDArray(0, i, xyz2, fragb[0][i])
            qcprot.SetDArray(1, i, xyz2, fragb[1][i])
            qcprot.SetDArray(2, i, xyz2, fragb[2][i])

        rot = qcprot.MakeDvector(9)
        rmsd = qcprot.CalcRMSDRotationalMatrix(xyz1, xyz2, fraglen, rot)
        if rmsd <= rmsd_cutoff:
            rotmat =[None] * 9
            for i in range(0,9):
                rotmat[i] = qcprot.GetDvector(i, rot)

            # translate the cluster to the smotif_coors:
            cluster_coo = formatClusterCoo(data)
            spin_coors = extractSpinCoo(cluster_coo, spin, res_type)
            cm_spin_coors = translateCM(spin_coors, b_cen)
            rot_spin_coors = applyRot(cm_spin_coors, rotmat)
            trans_spin_coors = applyTranslation(rot_spin_coors, a_cen)
            all_spin_coors = extend_array(all_spin_coors,  trans_spin_coors)
            #print trans_spin_coors
            #print all_spin_coors
            """
            cluster_coo = formatClusterCoo(data)
            cm_cluster_coo = translateCM(cluster_coo, b_cen)
            rot_cluster_coo = applyRot(cm_cluster_coo, rotmat)
            trans_cluster_coo = applyTranslation(rot_cluster_coo, a_cen)
            spin_coors = extractSpinCoo(trans_cluster_coo, spin, res_type )

            qcprot.FreeDMatrix(xyz1)
            qcprot.FreeDMatrix(xyz2)
            qcprot.FreeDArray(rot)

            #return trans_spin_coors
            """

        qcprot.FreeDMatrix(xyz1)
        qcprot.FreeDMatrix(xyz2)
        qcprot.FreeDArray(rot)

    if all_spin_coors:
        return all_spin_coors
    else:
        return False



