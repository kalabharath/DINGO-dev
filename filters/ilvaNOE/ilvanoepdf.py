from utility.io_util import readPickle
import bbRMSD

def getHCoorMatrix(ss_list, smotif):
    # change this to match the format from bbrmsd
    noe_matrix = {}
    count = 0
    for i in range(1, len(smotif), 5):
        coo1 = [smotif[i][3], smotif[i][4], smotif[i][5]]
        noe_matrix[ss_list[count]] = coo1
        count += 1
    return noe_matrix

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
    spin_coors = bbRMSD.bbrmsd(bbc, rotamers, rmsd_cutoff, spin, res_type )

    return spin_coors

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
    cluster_protons={}

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
        atom1_coor, atom2_coor = 0 , 0
        if (noedef[0] in resi) and (noedef[2] in resi):
            if noedef[1] == 'H':
                atom1_coor = coorH_matrix[noedef[0]]
            elif noedef[1] in methyls[noedef[7]]:
                if (noedef[0], noedef[1]) in cluster_proton_entries:
                    atom2_coor = cluster_protons[(noedef[0], noedef[1])]
                else:
                    bb1_coors = bb_matrix[noedef[0]]
                    atom2_coor = getILVARotamers(noedef[7], bb1_coors, noedef[1])
                    cluster_protons[(noedef[0], noedef[1])] = atom2_coor
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
                    cluster_protons[(noedef[2], noedef[3])] = atom2_coor
            else:
                print "Not a proton spin or missing residue number"
        else:
            pass

        if atom1_coor and atom2_coor:
            print noedef
            print atom1_coor, atom2_coor
            pass

    print "Done"
    return False
