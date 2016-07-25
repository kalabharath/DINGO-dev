import math
import numpy as np

"""
Module to provide subroutines to parse and fit RDC data
"""


def ZYZRot(A, B, G, scal=1.0):
    """
    *** Adopted from PyParaTools ***
    Returns a ZYZ rotation matrix when given 3 Euler angles (in degrees).
    See: http://mathworld.wolfram.com/EulerAngles.html
    :param A: The (A)lpha angle (degrees)
    :type  A: float
    :param B: The (B)eta angle (degrees)
    :type  B: float
    :param G: The (G)amma angle (degrees)
    :type  G: float
    :param scal: Scale the rotation matix by a constant [OPTIONAL]
    :type  scal: float
    :rtype: Numpy matrix
    """
    from math import cos
    from math import sin
    from math import radians

    rot = np.zeros((3, 3))
    ca = cos(radians(A))
    cb = cos(radians(B))
    cg = cos(radians(G))
    sa = sin(radians(A))
    sb = sin(radians(B))
    sg = sin(radians(G))
    rot[0][0] = ((-sg * sa) + (cb * ca * cg)) * scal
    rot[0][1] = ((sg * ca) + (cb * sa * cg)) * scal
    rot[0][2] = ((-cg * sb)) * scal
    rot[1][0] = ((-cg * sa) - (cb * ca * sg)) * scal
    rot[1][1] = ((cg * ca) - (cb * sa * sg)) * scal
    rot[1][2] = ((sg * sb)) * scal
    rot[2][0] = ((sb * ca)) * scal
    rot[2][1] = ((sb * sa)) * scal
    rot[2][2] = (cb) * scal
    return rot


def getGMR(spin_type1, spin_type2):
    """
    *** Adopted from PyParaTools ***
    Return the gyromagnetic ratio(s) for the given spin type.
    From: http://nmrwiki.org/wiki/index.php?title=Gyromagnetic_ratio
    The values defined are consistent with those in Xplor-NIH.
    :param spin_type: The atom identifier for the spin (can be H, C, N)
    :type spin_type: String (H, C, N)
    :rtype: float
    """
    from math import pi
    mgr = {'H': ((2 * pi) * 42.576) * 1e6, \
           'C': ((2 * pi) * 10.705) * 1e6, \
           'CA': ((2 * pi) * 10.705) * 1e6, \
           'N': ((2 * pi) * -4.315) * 1e6,}

    return mgr[spin_type1] * mgr[spin_type2]


def rdcScal(S, B0, temp):
    """
    Scaling constant.for RDC calculations
    """
    # TODO: These need to be checked
    hbar = 1.05457148e-34
    kboltz = 1.3806503e-23
    scal = -S * hbar * B0 * B0 / (8 * 15 * math.pi * math.pi * kboltz * temp)
    return scal * 0.01


def getVector(c1, c2):
    """
    return the vector for given two coordinates
    """

    vector = [c1[0][0] - c2[0][0], c1[0][1] - c2[0][1], c1[0][2] - c2[0][2]]
    return vector


def RDC_ZYZ(p0, scal, vec_data):
    """
    The RDC function for ZYZ notation
    :param p0:                     Parameters Dax, Drh, A, B, G
    :vec_data:                     Vector def, GMR*GMR, RDC
    :param scal:                   The RDC scaling constants
    """

    n = len(vec_data)
    err = np.zeros(n)
    rot = ZYZRot(p0[2], p0[3], p0[4])

    for i in xrange(n):
        # [[-0.6490001678466797, -0.3989999294281006, 0.6469993591308594], -7252794860688272.0, -6.593]]
        #             X                  Y                    Z                 gm1*gm2            rdc
        Vec = vec_data[i][0]
        GMRp = vec_data[i][1]
        rdc_measured = vec_data[i][2]
        kscal = scal * GMRp
        X = Vec[0]
        Y = Vec[1]
        Z = Vec[2]
        x_t = rot[0][0] * X + rot[0][1] * Y + rot[0][2] * Z
        y_t = rot[1][0] * X + rot[1][1] * Y + rot[1][2] * Z
        z_t = rot[2][0] * X + rot[2][1] * Y + rot[2][2] * Z
        r2 = (x_t * x_t) + (y_t * y_t) + (z_t * z_t)
        r5 = (r2 * r2) * math.sqrt(r2)
        tmp = 1.0 / r5
        rdc = ((p0[0] * (3.0 * z_t * z_t - r2) + p0[1] * 1.5 * (x_t * x_t - y_t * y_t)) * tmp) * kscal
        err[i] = rdc_measured - rdc
    return err


def RDCAxRhFit(s1_def, s2_def, smotif, exp_data):
    """

    :param s1_def:
    :param s2_def:
    :param smotif:
    :param threshold:
    :return:
    """

    ss1_list = range(s1_def[4], s1_def[5] + 1)
    ss2_list = range(s2_def[4], s2_def[5] + 1)




    return True