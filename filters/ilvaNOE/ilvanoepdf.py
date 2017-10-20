import bbRMSD
import math

def getHCoorMatrix(ss_list, smotif):
    # change this to match the format from bbrmsd
    noe_matrix = {}
    count = 0
    for i in range(1, len(smotif), 5):
        coo1 = [[smotif[i][3]], [smotif[i][4]], [smotif[i][5]], ['H']]
        noe_matrix[ss_list[count]] = coo1
        count += 1
    return noe_matrix

def processBBC(bbc):
    tlen = len(bbc)
    x, y, z = [None] * tlen, [None] * tlen, [None] * tlen
    for i in range (0, tlen):
        x[i] = bbc[i][3]
        y[i] = bbc[i][4]
        z[i] = bbc[i][5]
    return [x, y, z]

def getBackboneCoors(ss_list, smotif):
    backbone= {}
    count = 0
    j = 0
    for i in range(0, len(ss_list)):
        backbone[ss_list[count]] = smotif[j:j+5]
        j = j+5
        count += 1
    return backbone


def getILVARotamers(res_type, bbc, spin):
    import glob, os
    cwd = (os.path.dirname(os.path.realpath(__file__)))
    file_name = cwd+'/sidechainRotamers/'+res_type+"_sc/*.pickle"
    rotamers = glob.glob(file_name)
    rmsd_cutoff = 0.1
    bbc = processBBC(bbc)
    spin_coors = bbRMSD.bbrmsd(bbc, rotamers, rmsd_cutoff, spin, res_type )

    return spin_coors

def checkNoe(atom1_coor, atom2_coor, noedef):
    noe_bool = False
    dist = 999.999
    for i in range(0, len(atom1_coor[0])):
        x1, y1, z1 = atom1_coor[0][i], atom1_coor[1][i], atom1_coor[2][i]

        for j in range(0, len(atom2_coor[0])):
            x2, y2, z2 = atom2_coor[0][j], atom2_coor[1][j], atom2_coor[2][j]
            dist = math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1))
            if dist <= (noedef[4] + noedef[6]):
                noe_bool = True
                return noe_bool, dist
    return noe_bool, dist

def extractUnSatisfiedNoes(noes_found, data):
    unsatisfied = []
    for entry in data:
        if entry in noes_found:
            pass
        else:
            unsatisfied.append(entry)
    return unsatisfied
def s1ILVApdf(s1_def, s2_def, smotif, exp_data):
    """

    :param s1_def:
    :param s2_def:
    :param smotif:
    :param exp_data:
    :return:
    """

    methyls = {'I': ['HG2', 'HD1'], 'L': ['HD1', 'HD2', 'HD'], 'V': ['HG1', 'HG2', 'HG'], 'A': ['HB']}
    noe_data = exp_data['ilva_noes']
    cluster_protons = {}

    satisfied_noes = []
    unsatisfied_noes = []

    noes_found = 0.0
    total_noes = 0.0
    ss1_list = range(s1_def[4], s1_def[5] + 1)
    ss2_list = range(s2_def[4], s2_def[5] + 1)
    resi = ss1_list + ss2_list
    resi.sort()

    coorH_matrix = getHCoorMatrix(ss1_list, smotif[1])
    coor2_matrix = getHCoorMatrix(ss2_list, smotif[2])
    coorH_matrix.update(coor2_matrix)

    bb_matrix = getBackboneCoors(ss1_list, smotif[1])
    bb2_matrix = getBackboneCoors(ss2_list, smotif[2])
    bb_matrix.update(bb2_matrix)


    for noedef in noe_data:
        cluster_proton_entries = cluster_protons.keys()
        # noedef = [31, 'H', 47, 'HG1', 4.2, 2.4, 0.63, 'E', 'V']
        atom1_coor, atom2_coor = [] , []

        if (noedef[0] in resi) and (noedef[2] in resi):

            if noedef[1] == 'H':
                atom1_coor = coorH_matrix[noedef[0]]
            elif noedef[1] in methyls[noedef[7]]:
                if (noedef[0], noedef[1]) in cluster_proton_entries:
                    atom1_coor = cluster_protons[(noedef[0], noedef[1])]
                else:
                    bb1_coors = bb_matrix[noedef[0]]
                    atom1_coor = getILVARotamers(noedef[7], bb1_coors, noedef[1])
                    if atom1_coor:
                        cluster_protons[(noedef[0], noedef[1])] = atom1_coor
            else:
                print "Not a proton spin or missing residue number"

            if noedef[3] == 'H':
                atom2_coor = coorH_matrix[noedef[2]]
            elif noedef[3] in methyls[noedef[8]]:
                if (noedef[2], noedef[3]) in cluster_proton_entries:
                    atom2_coor = cluster_protons[(noedef[2], noedef[3])]
                else:
                    bb2_coors = bb_matrix[noedef[2]]
                    atom2_coor = getILVARotamers(noedef[8], bb2_coors, noedef[3])
                    if atom2_coor:
                        cluster_protons[(noedef[2], noedef[3])] = atom2_coor
            else:
                print "Not a proton spin or missing residue number"

            if atom1_coor and atom2_coor:
                noe_bool, dist = checkNoe(atom1_coor, atom2_coor, noedef)

                if noe_bool:

                    satisfied_noes.append(noedef)
                    noes_found += 1
                    total_noes += 1
                else:
                    total_noes += 1
            else:
                total_noes +=1
        else:
            pass
    unsatisfied_noes = extractUnSatisfiedNoes(satisfied_noes, noe_data)

    return (noes_found / total_noes), total_noes, [satisfied_noes, unsatisfied_noes], cluster_protons


def getSxCoorMatrix(coor_array, native_sse):

    resi ={}
    native_sse_range = range(native_sse[4], native_sse[5] + 1)
    count_array = 0
    for i in range(1, len(coor_array[0]), 5):
        resi[native_sse_range[count_array]] = [[coor_array[0][i]], [coor_array[1][i]], [coor_array[2][i]], ['H']]
        count_array += 1
    return resi


def getSxbbCoorMatrix(coor_array, native_sse):

    resi ={}
    native_sse_range = range(native_sse[4], native_sse[5] + 1)
    count_array = 0
    for i in range(0, len(coor_array[0]), 5):
        x, y, z = [None] * 5, [None] * 5, [None] * 5
        for j in range(0, 5):
            x[j] = coor_array[0][i + j]
            y[j] = coor_array[1][i + j]
            z[j] = coor_array[2][i + j]

        resi[native_sse_range[count_array]] = [x, y, z]
        count_array += 1
    return resi

def getSxILVARotamers(res_type, bbc, spin):
    import glob, os
    cwd = (os.path.dirname(os.path.realpath(__file__)))
    file_name = cwd+'/sidechainRotamers/'+res_type+"_sc/*.pickle"
    rotamers = glob.glob(file_name)
    rmsd_cutoff = 0.1

    spin_coors = bbRMSD.bbrmsd(bbc, rotamers, rmsd_cutoff, spin, res_type )

    return spin_coors
def sXILVApdf(transformed_coors, native_sse_order, current_ss, sorted_noe_data, cluster_protons):

    import copy
    sse_coors = copy.deepcopy(transformed_coors)
    satisfied_noes = sorted_noe_data[0][:]
    noe_data = sorted_noe_data[1]

    noes_found = len(satisfied_noes)
    total_noes = len(satisfied_noes)
    #print "noes_found, total_noes", noes_found, total_noes

    methyls = {'I': ['HG2', 'HD1'], 'L': ['HD1', 'HD2', 'HD'], 'V': ['HG1', 'HG2', 'HG'], 'A': ['HB']}
    coorH_matrix = {}
    #cluster_protons = {}

    for i in range(0, len(sse_coors)):
        tcoor = getSxCoorMatrix(sse_coors[i], native_sse_order[i])
        coorH_matrix.update(tcoor)
    resi = coorH_matrix.keys()

    bb_matrix = {}
    for i in range(0, len(sse_coors)):
        tcoor = getSxbbCoorMatrix(sse_coors[i], native_sse_order[i])
        bb_matrix.update(tcoor)

    for noedef in  sorted_noe_data[1] :

        cluster_proton_entries = cluster_protons.keys()
        # noedef = [31, 'H', 47, 'HG1', 4.2, 2.4, 0.63, 'E', 'V']
        atom1_coor, atom2_coor = [], []

        if (noedef[0] in resi) and (noedef[2] in resi):
            if noedef[1] == 'H':
                atom1_coor = coorH_matrix[noedef[0]]
            elif noedef[1] in methyls[noedef[7]]:
                if (noedef[0], noedef[1]) in cluster_proton_entries:
                    atom1_coor = cluster_protons[(noedef[0], noedef[1])]
                else:
                    bb1_coors = bb_matrix[noedef[0]]
                    atom1_coor = getSxILVARotamers(noedef[7], bb1_coors, noedef[1])
                    if atom1_coor:
                        cluster_protons[(noedef[0], noedef[1])] = atom1_coor
            else:
                print "Not a proton spin or missing residue number"

            if noedef[3] == 'H':
                atom2_coor = coorH_matrix[noedef[2]]
            elif noedef[3] in methyls[noedef[8]]:
                if (noedef[2], noedef[3]) in cluster_proton_entries:
                    atom2_coor = cluster_protons[(noedef[2], noedef[3])]
                else:
                    bb2_coors = bb_matrix[noedef[2]]
                    atom2_coor = getSxILVARotamers(noedef[8], bb2_coors, noedef[3])
                    if atom2_coor:
                        cluster_protons[(noedef[2], noedef[3])] = atom2_coor
            else:
                print "Not a proton spin or missing residue number"

            if atom1_coor and atom2_coor:
                noe_bool, dist = checkNoe(atom1_coor, atom2_coor, noedef)
                #print noe_bool, dist, noedef
                if noe_bool:
                    satisfied_noes.append(noedef)
                    noes_found += 1
                    total_noes += 1
                else:
                    total_noes += 1
            else:
                total_noes += 1
        else:
            pass

    #print "noes_found, total_noes", noes_found, total_noes
    unsatisfied_noes = extractUnSatisfiedNoes(satisfied_noes, noe_data)

    return (noes_found / total_noes), total_noes, [satisfied_noes, unsatisfied_noes], cluster_protons



