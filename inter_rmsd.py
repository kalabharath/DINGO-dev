import sys, os, copy
sys.path.append('../../main/')
__author__ = 'kalabharath'


import  utility.io_util as io

seq = int(sys.argv[1])


top_result = []

t_file = str(seq) + "_refined_tophits.gzip"
if os.path.isfile(t_file):
    top_result = io.readGzipPickle(t_file)

else:
    t_file = str(seq) + "_tophits.gzip"
    if os.path.isfile(t_file):
        top_result = io.readGzipPickle(t_file)
    else:
        print "Somethis is terrribly wrong !"
        exit()

for p in range(0, len(top_result)):
    print 'model_',p,

    top_struct = top_result[p]

    top_struct = copy.copy(top_struct)


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
    ss_list = top_struct[0][-1]
    exp_data = io.readPickle("exp_data.pickle")
    exp_data_types = exp_data.keys()  # ['ss_seq', 'pcs_data', 'aa_seq', 'contacts']

    aa_seq = exp_data['aa_seq']


