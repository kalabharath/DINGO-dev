import sys, os
sys.path.append("/home/kalabharath/zinr/main")
__author__ = 'kalabharath'


import  utility.io_util as io
import glob

def dumpPDBCoo2(coo_array, sse_pos, sse_seq, aa_seq,prefix):

    three_letter ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', \
                   'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR',    \
                   'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA',    \
                   'G':'GLY', 'P':'PRO', 'C':'CYS'}

    print sse_seq
    outfile = open(str(prefix)+"_"+str(sse_pos)+"_.pdb", 'w')
    t_count = 0
    res_num = sse_seq[-2]
    for i in range(0, len(coo_array[0])):

        x = coo_array[0][i]
        y = coo_array[1][i]
        z = coo_array[2][i]
        atom = coo_array[3][i]
        res_no = coo_array[4][i]
        res = coo_array[5][i]

        if t_count <= 5:
            t_count +=1
            if t_count < 5 :
                res = three_letter[aa_seq[res_num+1]]
                pdb_line = "%-6s%5d  %-2s%5s%2s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s"\
                   %('ATOM',i+1,atom,res,'A',res_num," ",x, y, z,1.0,30.0,' ',' \n')
                outfile.write(pdb_line)

            if t_count == 5:
                res = three_letter[aa_seq[res_num+1]]
                pdb_line = "%-6s%5d  %-2s%5s%2s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s"\
                   %('ATOM',i+1,atom,res,'A',res_num," ",x, y, z,1.0,30.0,' ',' \n')
                outfile.write(pdb_line)
                res_num+=1
                t_count = 0

        #pdb_line = "%-6s%5d  %-2s%5s%2s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s"\
        #           %('ATOM',i+1,atom,res,'A',res_no," ",x, y, z,1.0,30.0,' ',' \n')
        #outfile.write(pdb_line)
        #print pdb_line
    print i, sse_seq[-1]-sse_seq[-2]
    print 'TER'
    outfile.close()
    return True

def makeTopPickle3(previous_smotif_index, num_hits):

    hits = []
    regex=str(previous_smotif_index)+"_*_*.pickle"
    file_list = glob.glob(regex)
    for f in file_list:
        thits = io.readPickle(f)
        for thit in thits:
            hits.append(thit)
    """
                            0                                  1
    dump_log.append([transformed_coos,['seq_filter', csse_seq, seq_identity, blosum62_score],
                            2                                                           3
                ['contacts_filter', no_of_contacts, percent_of_satisfied_contacts], sse_ordered])

    """
    #print len(hits)
    new_dict={}
    seqs = []
    for hit in hits:
        #print hit
        seq_filter = hit[1]
        smotif_seq = seq_filter[1]
        contacts_filter = hit[2]
        if smotif_seq not in seqs:
            seqs.append(smotif_seq)
            dict_key_score = seq_filter[2]+contacts_filter[2]
            new_dict.setdefault(dict_key_score, []).append(hit)
    keys = new_dict.keys()
    keys.sort()
    dump_pickle = []
    try:
        for i in range(0,num_hits):
            dump_pickle.append(new_dict[keys[i]])
    except:
        print "could only make ", i
    io.dumpPickle(str(previous_smotif_index)+"_tophits.pickle", dump_pickle)
    return range(num_hits)


import utility.stage2_util as util
#makeTopPickle3(5, 10)

seq = int(sys.argv[1])


#util.makeTopPickle(seq,10)
import sys
top_result = io.readPickle(str(seq)+"_tophits.pickle")



for p in range(0, 10):

    top_struct = top_result[p]

    import copy

    top_struct = copy.copy(top_struct[0])

    print len(top_struct)

    sse_sequence =  top_struct[1][1]

    print sse_sequence

    coor_arrays = top_struct[2][1]
    #print coor_arrays
    ss_list = top_struct[0][-1]
    #print ss_list
    #print  len(coor_arrays)
    exp_data = io.readPickle("exp_data.pickle")
    exp_data_types = exp_data.keys()  # ['ss_seq', 'pcs_data', 'aa_seq', 'contacts']

    aa_seq = exp_data['aa_seq']

    print aa_seq
    for i in range(0, len(coor_arrays)):
        dumpPDBCoo2(coor_arrays[i], i, sse_sequence[i],aa_seq,p)


for q in range(0, 10):
    cat = "cat "+str(q)+"_*_.pdb >model_"+str(q)+".pdb"
    print cat
    os.system(cat)
