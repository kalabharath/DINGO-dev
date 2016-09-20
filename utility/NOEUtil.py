from numpy import zeros

import io_util as io


def getNOEData(noe_files, ss_seq):
    print noe_files
    noe_lines = io.readFile(noe_files)
    noe_matrix = zeros((len(ss_seq) + 1, len(ss_seq) + 1))
    for noel in noe_lines:

        if len(noel) <= 1:
            pass
        else:
            res1, atm1, res2, atom2, noe = noel.split()
            noe_matrix[int(res1), int(res2)] = float(noe)
            noe_matrix[int(res2), int(res1)] = float(noe)

    return noe_matrix
