#!/usr/bin/env python

"""
Project_Name: rmsd/qcp, File_name: qcp.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 8/05/15 , Time:4:27 AM
"""

import qcprot

def getCAcoo(frag):
    x,y,z = [], [], []
    for i in range(0, len(frag)):
        if frag[i][2] == 'CA':
            x.append(frag[i][3])
            y.append(frag[i][4])
            z.append(frag[i][5])
    return [x,y,z]

def rmsdQCP(psmotif, csmotif, direction):
    #print psmotif[0][0]
    #print csmotif[0][0]

    if direction =='left':
        frag_a = getCAcoo(psmotif[0][1])

    else:
        frag_a = getCAcoo(psmotif[0][2])

    if direction == 'left':
        frag_b = getCAcoo(csmotif[0][2])
    else:
        frag_b = getCAcoo(csmotif[0][1])

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

    """
     /* apply rotation matrix */
    for (int i = 0; i < len; ++i)
    {
        x = rotmat[0]*frag_b[0][i] + rotmat[1]*frag_b[1][i] + rotmat[2]*frag_b[2][i];
        y = rotmat[3]*frag_b[0][i] + rotmat[4]*frag_b[1][i] + rotmat[5]*frag_b[2][i];
        z = rotmat[6]*frag_b[0][i] + rotmat[7]*frag_b[1][i] + rotmat[8]*frag_b[2][i];

        frag_b[0][i] = x;
        frag_b[1][i] = y;
        frag_b[2][i] = z;
    }
    """
    for i in range(0,fraglen):
        x = rotmat[0]*frag_bt[0][i] + rotmat[1]*frag_bt[1][i] + rotmat[2]*frag_bt[2][i]
        y = rotmat[3]*frag_bt[0][i] + rotmat[4]*frag_bt[1][i] + rotmat[5]*frag_bt[2][i]
        z = rotmat[6]*frag_bt[0][i] + rotmat[7]*frag_bt[1][i] + rotmat[8]*frag_bt[2][i]
        frag_bt[0][i] = x
        frag_bt[1][i] = y
        frag_bt[2][i] = z

    print frag_bt

    print rmsd
    print frag_a
    return True, True
