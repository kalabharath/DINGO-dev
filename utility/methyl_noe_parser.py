def stripChar(atom):
    if "#" in atom:
        atom = atom.strip("#")
    if '*' in atom:
        atom = atom.strip("*")
    return atom


def checkMethyl(res, atom):
    methyls = {'I': ['HG2', 'HD1'], 'L': ['HD1', 'HD2', 'HD'], 'V': ['HG1', 'HG2', 'HG'], 'A': ['HB']}
    if "#" in atom:
        atom = atom.strip("#")
    if '*' in atom:
        atom = atom.strip("*")
    if res == 'A':
        if atom in methyls['A']:
            return True
        else:
            return False
    if res == 'V':
        if atom in methyls['V']:
            return True
        else:
            return False
    if res == 'I':
        if atom in methyls['I']:
            return True
        else:
            return False
    if res == 'L':
        if atom in methyls['L']:
            # print "different noes",res, atom,  methyls[res]
            return True
        else:
            return False

    return False


def identify_noe_type(tline, aa_seq, seq_correction):
    methyls = ['I', 'L', 'V', 'A']
    count_ors = 0
    for t in tline:
        if t == 'or':
            count_ors += 1

    if count_ors == 0:
        # ['assign', '(resid', '82', 'and', 'name', 'HN', ')', '(resid', '83', 'and', 'name', 'HN', ')', '2.5', '0.6', '0.6']
        #        print tline
        res1 = int(tline[2]) - seq_correction
        atom1 = tline[5].strip()
        res2 = int(tline[8]) - seq_correction
        atom2 = tline[11].strip()
        noe = float(tline[13])
        tol = float(tline[14])
        aa1 = aa_seq[res1 - 1]
        aa2 = aa_seq[res2 - 1]
        if res1 == res2:
            return []
        if abs(res2 - res1) < 3:
            return []
        if atom1 == 'HN' and atom2 == 'HN':
            return [res1, atom1, res2, atom2, noe, tol, aa1, aa2]
        elif aa1 in methyls and aa2 in methyls:
            atom1 = stripChar(atom1)
            check_aa1 = checkMethyl(aa1, atom1)
            if check_aa1:
                atom2 = stripChar(atom2)
                check_aa2 = checkMethyl(aa2, atom2)
                if check_aa2:
                    return [res1, atom1, res2, atom2, noe, tol, aa1, aa2]
        elif aa1 in methyls and atom2 == 'HN':
            atom1 = stripChar(atom1)
            check_aa1 = checkMethyl(aa1, atom1)
            if check_aa1:
                return [res1, atom1, res2, atom2, noe, tol, aa1, aa2]
        elif atom1 == 'HN' and aa2 in methyls:
            atom2 = stripChar(atom2)
            check_aa2 = checkMethyl(aa2, atom2)
            if check_aa2:
                return [res1, atom1, res2, atom2, noe, tol, aa1, aa2]
        else:
            return []
    elif count_ors == 2:
        # ['assign', '(resid', '84', 'and', '(name', 'HG1#', 'or', 'name', 'HG2#', '))', '(resid', '205', 'and', '(name', 'HD1#', 'or', 'name', 'HD2#', '))', '3.2', '0.8', '0.8']
        noe = float(tline[-3])
        tol = float(tline[-2])

        res1 = int(tline[2]) - seq_correction
        atom1 = tline[5].strip()
        atom12 = tline[8].strip()

        res2 = int(tline[11]) - seq_correction
        atom2 = tline[14].strip()
        atom22 = tline[17].strip()

        aa1 = aa_seq[res1 - 1]
        aa2 = aa_seq[res2 - 1]

        if res1 == res2:
            return []
        if abs(res2 - res1) < 3:
            return []
        if aa1 in methyls and aa2 in methyls:
            atom1 = stripChar(atom1)
            checkaa1 = checkMethyl(aa1, atom1)
            if checkaa1:
                atom12 = stripChar(atom12)
                checkatom12 = checkMethyl(aa1, atom12)
                if checkatom12:
                    atom2 = stripChar(atom2)
                    checkaa2 = checkMethyl(aa2, atom2)
                    if checkaa2:
                        atom22 = stripChar(atom22)
                        checkatom22 = checkMethyl(aa2, atom22)
                        if checkatom22:
                            return [[res1, atom1, res2, atom2, noe, tol, aa1, aa2],
                                    [res1, atom12, res2, atom2, noe, tol, aa1, aa2],
                                    [res1, atom1, res2, atom22, noe, tol, aa1, aa2],
                                    [res1, atom12, res2, atom22, noe, tol, aa1, aa2]]


    elif count_ors == 1:
        # ['assign', '(resid', '90', 'and', '(name', 'HB1', 'or', 'name', 'HB2', '))', '(resid', '92', 'and', 'name', 'HN', ')', '3.6', '0.9', '0.9']
        # ['assign', '(resid', '91', 'and', 'name', 'HD#', ')', '(resid', '94', 'and', '(name', 'HD2#', 'or', 'name', 'HD1#', '))', '2.7', '0.7', '0.7']

        noe = float(tline[-3])
        tol = float(tline[-2])

        if tline[6] == 'or':
            res1 = int(tline[2]) - seq_correction
            atom1 = tline[5].strip()
            atom12 = tline[8].strip()

            res2 = int(tline[11]) - seq_correction
            atom2 = tline[14].strip()
            aa1 = aa_seq[res1 - 1]
            aa2 = aa_seq[res2 - 1]
            if res1 == res2:
                return []
            if abs(res2 - res1) < 3:
                return []
            if aa1 in methyls and aa2 in methyls:
                atom1 = stripChar(atom1)
                checkaa1 = checkMethyl(aa1, atom1)
                if checkaa1:
                    atom12 = stripChar(atom12)
                    checkatom12 = checkMethyl(aa1, atom12)
                    if checkatom12:
                        atom2 = stripChar(atom2)
                        checkaa2 = checkMethyl(aa2, atom2)
                        if checkaa2:
                            return [[res1, atom1, res2, atom2, noe, tol, aa1, aa2],
                                    [res1, atom12, res2, atom2, noe, tol, aa1, aa2]]
            elif aa1 in methyls and atom2 == 'HN':
                atom1 = stripChar(atom1)
                checkaa1 = checkMethyl(aa1, atom1)
                if checkaa1:
                    atom12 = stripChar(atom12)
                    checkatom12 = checkMethyl(aa1, atom12)
                    if checkatom12:
                        return [[res1, atom1, res2, atom2, noe, tol, aa1, aa2],
                                [res1, atom12, res2, atom2, noe, tol, aa1, aa2]]
            else:
                return []

        if tline[12] == 'or':
            res1 = int(tline[2]) - seq_correction
            atom1 = tline[5].strip()

            res2 = int(tline[8]) - seq_correction
            atom2 = tline[11].strip()
            atom22 = tline[14].strip()

            aa1 = aa_seq[res1 - 1]
            aa2 = aa_seq[res2 - 1]
            if res1 == res2:
                return []
            if abs(res2 - res1) < 3:
                return []
            if aa1 in methyls and aa2 in methyls:
                atom1 = stripChar(atom1)
                checkaa1 = checkMethyl(aa1, atom1)
                if checkaa1:
                    atom2 = stripChar(atom2)
                    checkatom2 = checkMethyl(aa2, atom22)
                    if checkatom2:
                        atom22 = stripChar(atom22)
                        checkaa2 = checkMethyl(aa2, atom22)
                        if checkaa2:
                            return [[res1, atom1, res2, atom2, noe, tol, aa1, aa2],
                                    [res1, atom1, res2, atom22, noe, tol, aa1, aa2]]
            elif atom1 == 'HN' and aa2 in methyls:
                atom2 = stripChar(atom2)
                checkatom2 = checkMethyl(aa2, atom22)
                if checkatom2:
                    atom22 = stripChar(atom22)
                    checkaa2 = checkMethyl(aa2, atom22)
                    if checkaa2:
                        return [[res1, atom1, res2, atom2, noe, tol, aa1, aa2],
                                [res1, atom1, res2, atom22, noe, tol, aa1, aa2]]
            else:
                return []
    else:
        return []
    return []


def parseMethylNoes(seq_file, cyana_file, seq_correction):
    with open(seq_file) as fin:
        aa_lines = fin.readlines()
    aa_seq = ''
    for i in range(1, len(aa_lines)):
        taa = aa_lines[i].strip()
        aa_seq = aa_seq + taa

    with open(cyana_file) as fin:
        lines = fin.readlines()
    seq_correction = 78

    all_noes = []
    for line in lines:
        tline = line.split()
        if len(tline) > 1 and (tline[0].strip() == 'assign'):
            noes = identify_noe_type(tline, aa_seq, seq_correction)
            if noes:
                all_noes.append(noes)
    return all_noes

def methylNOEParser(block_file, aa_seq, seq_correction):
    print  aa_seq
    with open(block_file) as fin:
        lines = fin.readlines()
    all_noes = []
    for line in lines:
        tline = line.split()
        if len(tline) > 1 and (tline[0].strip() == 'assign'):
            noes = identify_noe_type(tline, aa_seq, seq_correction)
            if noes:
                all_noes.append(noes)
    return all_noes, len(all_noes)

def main():
    with open('seq.fasta') as fin:
        aa_lines = fin.readlines()
    aa_seq = ''
    for i in range(1, len(aa_lines)):
        taa = aa_lines[i].strip()
        aa_seq = aa_seq + taa

    with open('block_515450.mr') as fin:
        lines = fin.readlines()
    seq_correction = 78

    all_noes = []
    for line in lines:
        tline = line.split()
        if len(tline) > 1 and (tline[0].strip() == 'assign'):
            noes = identify_noe_type(tline, aa_seq, seq_correction)
            if noes:
                all_noes.append(noes)
    for n in all_noes:
        print n

    print len(all_noes)
    return True

if __name__ == '__main__':
    main()
