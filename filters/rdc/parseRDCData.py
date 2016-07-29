def FormatRdc(seqlen, rdcfile):
    """
    parses rdc from .npc file.
    should be in the format #['179', 'H', '179', 'N', '16.042', '0.0']
    the rdcs are returned as a dict with res_no as key and rdc def as value.
    """
    with open(rdcfile) as f:
        rdc_l = f.read().splitlines()
    rdcs = {}
    for l in rdc_l:
        r1, v1, r2, v2, rdc, tol = l.split()
        rdcs.setdefault(int(r1), []).append([int(r1), v1, int(r2), v2, float(rdc)])
    return rdcs

def getRdcData(rdc_in_file):
        """
        Parse RDCs from files
        the file names should be given inside a text file
        preceeded with flag '-in:file:rdc' as used in
        input for minirosetta.
        """
        with open(rdc_in_file,'r') as rdcin:
            files = rdcin.readlines()

        rdc_files = []
        for f in files:
            if(f[0:12] == '-in:file:rdc'):
                 farray = f.split()
                 for i in range(1, len(farray)):
                     rdc_files.append(farray[i])
        rdc_data = []
        for j in range(0, len(rdc_files)):
            rdc_data.append(FormatRdc(len(self.getDataSeq()), rdc_files[j]))
        return rdc_data
