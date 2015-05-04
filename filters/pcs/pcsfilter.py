#!/usr/bin/env python

"""
Project_Name: main/filters/pcsfilter, File_name: pcsfilter.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 24/04/15 , Time: 01:30 PM

fit PCS and apply Ax,Rh filters
"""

import fastT1FM

def getHN(ss1_list, ss2_list, smotif, atom_type='H'):
    """
    parse the x,y,z of individual SSE into seperate arrays
    :param smotif:
    :return: arrays of coordinates of individual SSE in seperate arrays
    """
    #TODO return arrays for any given atom type(s)

    rH1 , rH2 = [], []
    counter = 0
    for entry in  smotif[0][1]:
        # res_no, aa_type, atom_type, x, y, z
        if entry[2] == atom_type:
            x, y, z = entry[3], entry[4], entry[5]
            res_no = ss1_list[counter]
            counter +=1
            if counter == len(ss1_list):
                break
            rH1.append([x, y, z,res_no])
    counter = 0
    for entry in  smotif[0][2]:
        # res_no, aa_type, atom_type, x, y, z
        if entry[2] == atom_type:
            x, y, z = entry[3], entry[4], entry[5]
            res_no = ss2_list[counter]
            counter +=1
            if counter == len(ss2_list):
                break
            rH1.append([x, y, z,res_no])
    return rH1, rH2


def match_pcss_HN(rh1, rh2, pcs_data):
    """

    :param rh1:
    :param rh2:
    :param pcs_data:
    :return:
    """
    smotif_pcs = []

    for entry in rh1:
        res_no = (entry[-1])-1
        smotif_pcs.append(pcs_data[res_no])
    for entry in rh2:
        res_no = (entry[-1])-1
        smotif_pcs.append(pcs_data[res_no])
    return  smotif_pcs

def PointsOnSpheres(M, N, rMx, rMy, rMz):
    """
    quick way from wikipedia
    :param M:
    :param N:
    :param rMx:
    :param rMy:
    :param rMz:
    :return:
    """

    import math
    node = []
    dlong = math.pi*(3-math.sqrt(5))  # ~2.39996323
    dz = 2.0/N
    xlong = 0
    z = 1 - 0.5*dz
    for k in range(N):
        r = math.sqrt(1-z*z)
        node.append([math.cos(xlong)*r, math.sin(xlong)*r, z])
        z -= dz
        xlong += dlong

    j = 0
    for i in range(M[0],M[1]):
        for k in range(N):
            fastT1FM.SetDvector(j, rMx, i*node[k][0])
            fastT1FM.SetDvector(j, rMy, i*node[k][1])
            fastT1FM.SetDvector(j, rMz, i*node[k][2])
            j += 1

def usuablePCS(pcs_array):

    for j in range(0, len(pcs_array[0])):
        counter = 0
        for entry in pcs_array:
            if entry[j] != 999.999:
                counter +=1
        if counter <=5:
            return False
    return True

def calc_axrh(saupe_matrices):
    """

    :param saupe_matrices:
    :return:
    """
    import numpy
    from numpy import linalg as LA
    axrh=[]
    for t in saupe_matrices:
        w, v = LA.eig(numpy.array([[t[0], t[1], t[2]], [t[1], t[3], t[4]], [t[2], t[4], -t[0]-t[3]]]))
        x = []
        for i in range(3):
            x.append([abs(w[i]), w[i]])
        x.sort()
        for i in range(3):
            w[i] = x[i][1]

        axrh.append([w[2]-0.5*(w[0]+w[1]),w[0]-w[1]])
    return axrh

@profile
def PCSAxRhFit(s1_def, s2_def, smotif, exp_data, threshold = 0.05):
    """

    :param s1_def:
    :param s2_def:
    :param smotif:
    :param threshold:
    :return:
    """

    ss1_list = range(s1_def[4], s1_def[5]+1)
    ss2_list = range(s2_def[4], s2_def[5]+1)

    smotif_ss1 = range(int(smotif[0][0][1]), int(smotif[0][0][2])+1 )
    smotif_ss2 = range(int(smotif[0][0][3]), int(smotif[0][0][4])+1 )


    #print ss1_list, ss2_list
    #print smotif_ss1, smotif_ss2
    #print smotif[0][0]

    rH1, rH2 = getHN(ss1_list, ss2_list, smotif, atom_type='H')
    pcs_data = exp_data['pcs_data']
    ntags = len(pcs_data)

    # Define Thomas's implementaion of hollow concentric shells

    nM = 1000  # 1000 pts in each sphere
    M = [1,40] # 40 spheres 10-50 Angstrom
    npts = (M[1] - M[0]) * nM  # 50 spheres * 1000 pts each
    rMx = fastT1FM.MakeDvector(npts) #allocate memmory
    rMy = fastT1FM.MakeDvector(npts)
    rMz = fastT1FM.MakeDvector(npts)
    PointsOnSpheres(M, nM, rMx, rMy, rMz)


    for tag in range(0,ntags):
    #for tag in range(0,1):
        smotif_pcs = match_pcss_HN(rH1, rH2, pcs_data[tag])

        # initialize variables for Thomas's PCS fitting routine
        frag_len = len(smotif_pcs)
        nsets = len(smotif_pcs[0])
        xyz = fastT1FM.MakeDMatrix(frag_len, 3)
        pcs = fastT1FM.MakeDMatrix(nsets, frag_len)
        xyz_HN = rH1+rH2

        if len(xyz_HN) != frag_len:
            print xyz_HN, len(xyz_HN)
            print smotif_pcs, len(smotif_pcs)


        for k in range(nsets):
            for j in range(frag_len):
                fastT1FM.SetDArray(k, j, pcs, smotif_pcs[j][k])

        cm = [0.0, 0.0, 0.0]
        for j in range(frag_len):
            cm[0] = cm[0] + xyz_HN[j][0]
            cm[1] = cm[1] + xyz_HN[j][1]
            cm[2] = cm[2] + xyz_HN[j][2]
        cm[0] /= float(frag_len)
        cm[1] /= float(frag_len)
        cm[2] /= float(frag_len)
        for j in range(frag_len):
            fastT1FM.SetDArray(j, 0, xyz, xyz_HN[j][0]-cm[0])
            fastT1FM.SetDArray(j, 1, xyz, xyz_HN[j][1]-cm[1])
            fastT1FM.SetDArray(j, 2, xyz, xyz_HN[j][2]-cm[2])

        tensor = fastT1FM.MakeDMatrix(nsets, 8)
        #ttmp = fastT1FM.MakeDvector(5)
        #Xtmp = fastT1FM.MakeDvector(2)
        Xaxrh_range = fastT1FM.MakeDMatrix(nsets, 4)
        for i in range(0,nsets):
            fastT1FM.SetDArray(i, 0, Xaxrh_range, 0.05)
            fastT1FM.SetDArray(i, 1, Xaxrh_range, 100.0)
            fastT1FM.SetDArray(i, 2, Xaxrh_range, 0.05)
            fastT1FM.SetDArray(i, 3, Xaxrh_range, 100.0)

        if usuablePCS(smotif_pcs):

            chisqr = fastT1FM.rfastT1FM_multi(npts, rMx, rMy, rMz, nsets, frag_len, xyz, pcs, tensor, Xaxrh_range)

        else:
            chisqr = 1.0e+30

        if chisqr < 1.0e+30:
            x = fastT1FM.GetDArray(0, 0, tensor)
            y = fastT1FM.GetDArray(0, 1, tensor)
            z = fastT1FM.GetDArray(0, 2, tensor)
            saupe_array=[]
            for kk in range(nsets):
                temp_saupe=[]
                for j in range(3,8):
                    temp_saupe.append(fastT1FM.GetDArray(kk, j, tensor))
                    #fastT1FM.SetDvector(j-3, ttmp, fastT1FM.GetDArray(kk, j, tensor))
                saupe_array.append(temp_saupe)
            metalpos=[x+cm[0], y+cm[1], z+cm[2]]
            AxRh = calc_axrh(saupe_array)
            #print tag+1, chisqr, metalpos, AxRh
    return True