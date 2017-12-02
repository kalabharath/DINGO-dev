import sys, os
#sys.path.append('/short/xc4/kbp502/BOSS/zinr')
sys.path.append('../../main/')
__author__ = 'kalabharath'


import  utility.io_util as io
import glob


import utility.stage2_util as util


seq = int(sys.argv[1])
num_hits = 5000
stage = 4

top_result = io.readGzipPickle(str(seq)+"_tophits.gzip")


for p in range(0, len(top_result)):
#for p in range(0,5):
    print 'model_',p,
    top_struct = top_result[p]

    import copy

    top_struct = copy.copy(top_struct)
    #print len(top_struct)

    sse_sequence =  top_struct[1][1]

    #print sse_sequence
    for entry in top_struct:
        if entry[0] =='cathcodes':
            print entry
            pass
        if entry[0] == 'Ref_RMSD':
            print entry[:-1]
        if entry[0] == 'RDC_filter':
            pass
            print entry
        if entry[0] == 'NOE_filter':
            print entry[0:4]
    #coor_arrays = top_struct[2][1]
    #print top_struct[2][0]
    #print "no of sse elements", len(coor_arrays)
    ss_list = top_struct[0][-1]
    #print ss_list
    #print  len(coor_arrays)
    exp_data = io.readPickle("exp_data.pickle")
    exp_data_types = exp_data.keys()  # ['ss_seq', 'pcs_data', 'aa_seq', 'contacts']

    aa_seq = exp_data['aa_seq']

    #print aa_seq, len(aa_seq)
