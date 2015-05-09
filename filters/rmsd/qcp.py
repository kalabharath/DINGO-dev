#!/usr/bin/env python

"""
Project_Name: rmsd/qcp, File_name: qcp.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 8/05/15 , Time:4:27 AM
"""

import qcprot

def dumpPDBCoo(coo_array, atom_type='CA'):
    for i in range(0, len(coo_array[0])):
        x = coo_array[0][i]
        y = coo_array[1][i]
        z = coo_array[2][i]
        pdb_line =



    """
    tline="%-6s%5d %2s%1s%s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s"%\
          ('ATOM',i,line[12:16],"",res_name,'A',res_no+1," ",float(line[30:38]),float(line[38:46]),float(line[46:54]),1\
          .0,30.0,' ',' ')
    """
    pass
    return True

def getCAcoo(frag):
    x,y,z = [], [], []
    for i in range(0, len(frag)):
        if frag[i][2] == 'CA':
            x.append(frag[i][3])
            y.append(frag[i][4])
            z.append(frag[i][5])
    return [x,y,z]

def getcoo(frag):
    x,y,z = [], [], []
    for i in range(0, len(frag)):
        x.append(frag[i][3])
        y.append(frag[i][4])
        z.append(frag[i][5])
    return [x,y,z]

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


def rmsdQCP(psmotif, csmotif, direction):
    #print psmotif[0][0]
    #print csmotif[0][0]

    if direction =='left':
        native_frag_a = getCAcoo(psmotif[0][1])

    else:
        native_frag_a = getCAcoo(psmotif[0][2])

    if direction == 'left':
        native_frag_b = getCAcoo(csmotif[0][2])
    else:
        native_frag_b = getCAcoo(csmotif[0][1])

    frag_a, a_cen = centerCoo(native_frag_a)
    frag_b, b_cen = centerCoo(native_frag_b)

    fraglen = len(frag_a[0])

    xyz1 = qcprot.MakeDMatrix(3, fraglen)
    xyz2 = qcprot.MakeDMatrix(3, fraglen)

    for i in range(0, fraglen):
        qcprot.SetDArray(0, i, xyz1, frag_a[0][i])
        qcprot.SetDArray(1, i, xyz1, frag_a[1][i])
        qcprot.SetDArray(2, i, xyz1, frag_a[2][i])
    for i in range(0, fraglen):
        qcprot.SetDArray(0, i, xyz2, frag_b[0][i])
        qcprot.SetDArray(1, i, xyz2, frag_b[1][i])
        qcprot.SetDArray(2, i, xyz2, frag_b[2][i])

    rot = qcprot.MakeDvector(9)

    #*********
    rmsd = qcprot.CalcRMSDRotationalMatrix(xyz1, xyz2, fraglen, rot)
    #*********

    tx, ty, tz = [], [], []
    for i in range(0, fraglen):
        tx.append(qcprot.GetDArray(0, i, xyz2))
        ty.append(qcprot.GetDArray(1, i, xyz2))
        tz.append(qcprot.GetDArray(2, i, xyz2))

    frag_bt = [tx, ty, tz]
    print frag_bt

    rotmat = []
    for i in range(0,9):
        rotmat.append(qcprot.GetDvector(i,rot))

    for i in range(0,fraglen):
        x = rotmat[0]*frag_bt[0][i] + rotmat[1]*frag_bt[1][i] + rotmat[2]*frag_bt[2][i]
        y = rotmat[3]*frag_bt[0][i] + rotmat[4]*frag_bt[1][i] + rotmat[5]*frag_bt[2][i]
        z = rotmat[6]*frag_bt[0][i] + rotmat[7]*frag_bt[1][i] + rotmat[8]*frag_bt[2][i]
        frag_bt[0][i] = x
        frag_bt[1][i] = y
        frag_bt[2][i] = z

    print rmsd

    return True, True
