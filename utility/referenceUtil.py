import utility.io_util as io


def getRefCoors(pdb_file):
    """
    search for 'CAs'
    :param pdb_file:
    :return:
    """

    ca_coors = []
    data = io.readFile(pdb_file)
    for line in data:

        if line[12:16].strip() == 'CA':
            xcoor = float(line[30:38])
            ycoor = float(line[38:46])
            zcoor = float(line[46:54])
            ca_coors.append([xcoor, ycoor, zcoor])

    return ca_coors