#!/usr/bin/env python

"""
Project_Name: rmsd/qcp, File_name: qcp.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 8/05/15 , Time:4:27 AM
"""

import qcprot

def dumpPDBCoo(coo_array):
    for i in range(0, len(coo_array[0])):
        x = coo_array[0][i]
        y = coo_array[1][i]
        z = coo_array[2][i]

        pdb_line = "%-6s%5d %2s%1s%s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s"\
                   %('ATOM',i+1,'CA',"",'ALA','A',i+1," ",x, y, z,1.0,30.0,' ',' ')
        print pdb_line
    print 'TER'
    return True

def dumpPDBCoo2(coo_array):
    for i in range(0, len(coo_array[0])):
        x = coo_array[0][i]
        y = coo_array[1][i]
        z = coo_array[2][i]
        atom = coo_array[3][i]
        res_no = coo_array[4][i]
        res = coo_array[5][i]

        pdb_line = "%-6s%5d  %-2s%5s%2s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s"\
                   %('ATOM',i+1,atom,res,'A',res_no," ",x, y, z,1.0,30.0,' ',' ')
        print pdb_line
    print 'TER'
    return True

def getCAcoo(frag):
    #print frag
    x,y,z = [], [], []
    for i in range(0, len(frag[0])):
        if frag[3][i]=='CA':
            x.append(frag[0][i])
            y.append(frag[1][i])
            z.append(frag[2][i])
    return [x,y,z]

def getcoo(frag):
    #[41, 'ASP', 'N', 28.117, -17.694, 7.215]
    x,y,z, atom_type, res_no, res = [], [], [], [], [], []
    for i in range(0, len(frag)):
        x.append(frag[i][3])
        y.append(frag[i][4])
        z.append(frag[i][5])
        atom_type.append(frag[i][2])
        res_no.append(frag[i][0])
        res.append(frag[i][1])
    return [x,y,z, atom_type, res_no, res]

def centerCoo(coo_array):

    xsum, ysum, zsum = 0,0,0

    for i in range(0, len(coo_array[0])):
        xsum += coo_array[0][i]
        ysum += coo_array[1][i]
        zsum += coo_array[2][i]

    xsum /=len(coo_array[0])
    ysum /=len(coo_array[0])
    zsum /=len(coo_array[0])

    for i in range(0, len(coo_array[0])):
        coo_array[0][i] -= xsum
        coo_array[1][i] -= ysum
        coo_array[2][i] -= zsum
    return coo_array, [xsum, ysum, zsum]

def applyTranslation(frag, cen_mass):

    for i in range(0, len(frag[0])):
        frag[0][i] += cen_mass[0]
        frag[1][i] += cen_mass[1]
        frag[2][i] += cen_mass[2]
    return frag

def translateCM(coo_array, cm):
    xsum, ysum, zsum = cm[0], cm[1], cm[2]
    for i in range(0, len(coo_array[0])):
        coo_array[0][i] -= xsum
        coo_array[1][i] -= ysum
        coo_array[2][i] -= zsum
    return coo_array

def applyRot(frag, rotmat):
    for i in range(0,len(frag[0])):
        x = rotmat[0]*frag[0][i] + rotmat[1]*frag[1][i] + rotmat[2]*frag[2][i]
        y = rotmat[3]*frag[0][i] + rotmat[4]*frag[1][i] + rotmat[5]*frag[2][i]
        z = rotmat[6]*frag[0][i] + rotmat[7]*frag[1][i] + rotmat[8]*frag[2][i]
        frag[0][i] = x
        frag[1][i] = y
        frag[2][i] = z
    return frag


def rmsdQCP(psmotif, csmotif, direction):
    #print psmotif[0][0]
    #print csmotif[0][0]

    if direction =='left':
        native_frag_a = psmotif[0][1]
        frag_a = getcoo(native_frag_a)
        native_frag_b = csmotif[0][2]
        frag_b = getcoo(native_frag_b)
    else:

        native_frag_a = psmotif[0][2]
        frag_a = getcoo(native_frag_a)
        native_frag_b = csmotif[0][1]
        frag_b = getcoo(native_frag_b)

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

    #*********
    rmsd = qcprot.CalcRMSDRotationalMatrix(xyz1, xyz2, fraglen, rot)
    #*********

    """
    tx, ty, tz = [], [], []
    for i in range(0, fraglen):
        tx.append(qcprot.GetDArray(0, i, xyz2))
        ty.append(qcprot.GetDArray(1, i, xyz2))
        tz.append(qcprot.GetDArray(2, i, xyz2))

    frag_bt = [tx, ty, tz]
    #print frag_bt
    """
    rotmat = []
    for i in range(0,9):
        rotmat.append(qcprot.GetDvector(i,rot))

    #dumpPDBCoo(frag_bca)

    rotated_fragb = applyRot(frag_b, rotmat)

    #dumpPDBCoo2(rotated_fragb)
    trans_fragb = applyTranslation(rotated_fragb, a_cen)
    #dumpPDBCoo2(trans_fragb)

    #print rmsd

    # translate the other SSE of the current smotif

    if direction == 'left':
        sse_2nd = csmotif[0][1]
    else:
        sse_2nd = csmotif[0][2]
    #print sse_2nd

    sse_2nd_coos = getcoo(sse_2nd)
    cm_sse2nd = translateCM(sse_2nd_coos,b_cen)
    rot_sse_2nd = applyRot(cm_sse2nd, rotmat)
    trans_sse2nd = applyTranslation(rot_sse_2nd, a_cen)
    dumpPDBCoo2(trans_sse2nd)

    return True, True
