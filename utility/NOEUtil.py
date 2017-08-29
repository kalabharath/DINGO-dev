from numpy import zeros

import io_util as io
import methyl_noe_parser

def getNOEData(noe_files, ss_seq):
    """
    Depreciated
    :param noe_files:
    :param ss_seq:
    :return:
    """
    print noe_files

    noe_lines = io.readFile(noe_files)
    noe_matrix = zeros((len(ss_seq) + 1, len(ss_seq) + 1))
    total_noe_count = 0
    for noel in noe_lines:

        if len(noel) <= 1:
            pass
        else:
            res1, atm1, res2, atom2, noe = noel.split()
            total_noe_count += 1
            noe_matrix[int(res1), int(res2)] = float(noe)
            noe_matrix[int(res2), int(res1)] = float(noe)

    return noe_matrix, total_noe_count


def parseNOEData(noe_files):
    """
    New NOE scoring function
    :param noe_files:
    :return:
    """
    print noe_files
    noe_lines = io.readFile(noe_files)
    noe_data = []
    total_noe_count = 0.0
    for noe in noe_lines:
        if len(noe) <= 1:
            pass
        else:
            res1, atm1, res2, atom2, noe, tol = noe.split()
            total_noe_count += 1.0
            if float(tol) == 0.0:
                tol = 0.1
            noe_data.append([int(res1), atm1, int(res2), atom2, float(noe), float(tol)])
    return noe_data, total_noe_count

def parseBMRBblockMR(block_file, aa_seq, seq_correction):

    return methyl_noe_parser.methylNOEParser(block_file, aa_seq, seq_correction)
