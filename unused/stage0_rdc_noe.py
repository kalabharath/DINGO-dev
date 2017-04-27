#!/usr/bin/env python

"""
Project_Name: main, File_name: stage0_rdc_noe.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 08/08/2016

prepare data to run with BOSS-R
"""

import sys

sys.path.append("../../main")
import multiprocessing

import utility.RDCUtil    as ru
import utility.io_util    as io
import utility.ss_util    as ss
import utility.NOEUtil    as nu
import utility.referenceUtil as ref

# Parse the input data text file
data = io.readInputDataFiles('input_data.txt')
datatypes = data.keys()
print datatypes

# Read in the amino acid sequence and the Secondary structure assignment
handle, aa_seq = io.readFasta(data['fasta_file'])
ss_seq = ru.matchSeq2SS(aa_seq, data['ss_file'])
print ss_seq

# Generate fuzzy +/-2 SSE combinations
ss_def, ss_combi = ss.genSSCombinations(ss_seq)
io.dumpPickle("ss_profiles.pickle", ss_combi)
print ss_def
# Read the native pdbs that you can exclude from the smotif search
native_pdbs = data['native_pdbs']
native_pdbs = native_pdbs.lower()
native_pdbs = native_pdbs.split()
print native_pdbs

try:

    homologs = data['homolog_pdbs']
    homologs = homologs.lower()
    homologs = homologs.split()
    thomologs = []
    for entry in homologs:
        if len(entry) > 4:
            thomologs.append(entry[0:4])
        else:
            thomologs.append(entry)
    print thomologs
except:
    pass

pred_axial = data['predicted_axial']
pred_axial = pred_axial.split()
pred_axial = [float(i) for i in pred_axial]

TinK = data['TinK']
TinK = TinK.split()
TinK = [float(i) for i in TinK]
print TinK

B0inT = data['B0inT']
B0inT = B0inT.split()
B0inT = [float(i) for i in B0inT]
print B0inT

exp_error = data['exp_error']
exp_error = exp_error.split()
exp_error = [float(i) for i in exp_error]

abs_exp_error = data['abs_exp_error']
abs_exp_error = abs_exp_error.split()
abs_exp_error = [float(i) for i in abs_exp_error]

rank_top_hits = data['rank_top_hits']
rank_top_hits = rank_top_hits.split()
rank_top_hits = [float(i) for i in rank_top_hits]

noe_fmeasure = data['noe_fmeasure']
noe_fmeasure = noe_fmeasure.split()
noe_fmeasure = [float(i) for i in noe_fmeasure]

rmsd_cutoff = data['rmsd_cutoff']
rmsd_cutoff = rmsd_cutoff.split()
rmsd_cutoff = [float(i) for i in rmsd_cutoff]

reference_ca = ref.getRefCoors(data['reference_pdb'])
print rmsd_cutoff

clash_distance = float(data['clash_distance'])
print 'clash_distance: ', clash_distance

rdc_data = ru.getRdcData(data['rdc_input_files'], ss_seq)
noe_data = nu.getNOEData(data['noe_input_files'], ss_seq)

map_route = ru.getRDCMapRoute(ss_combi, rdc_data)
print map_route

io.dumpPickle("rdc_route.pickle", map_route)

database_cutoff = data['database_cutoff']

data_keys = ['ss_seq', 'rdc_data', 'aa_seq', 'natives', 'clash_distance', 'database_cutoff', 'rmsd_cutoff',
             'reference_ca', 'pred_axial', 'exp_error', 'abs_exp_error', 'noe_data', 'noe_fmeasure', 'rank_top_hits'
                                                                                                     'B0inT', 'TinK',
             'homologs']

data_dict = {'ss_seq': ss_seq, 'rdc_data': rdc_data, 'aa_seq': aa_seq, 'natives': native_pdbs,
             'clash_distance': clash_distance, 'database_cutoff': database_cutoff,
             'rmsd_cutoff': rmsd_cutoff, 'reference_ca': reference_ca,
             'pred_axial': pred_axial, 'exp_error': exp_error, 'abs_exp_error': abs_exp_error, 'noe_data': noe_data,
             'noe_fmeasure': noe_fmeasure, 'rank_top_hits': rank_top_hits, 'B0inT': B0inT, 'TinK': TinK,
             }

io.dumpPickle("exp_data.pickle", data_dict)

fout = open("run.sh", 'w')
fout.write("#!/bin/bash\n")
ncpus = multiprocessing.cpu_count()

for i in range(0, len(map_route)):
    smotif = map_route[i]
    print smotif
    if i == 0:
        run_line = "mpirun -np " + str(ncpus) + " python ../../main/stage1_mpi_run.py\n"
        print run_line
        fout.write(run_line)
    elif i == 1:
        run_line = "mpirun -np " + str(ncpus) + " python ../../main/stage2_mpi_run.py 1000\n"
        print run_line
        fout.write(run_line)
    elif i != 1 and i <= len(map_route) - 3:
        run_line = "mpirun -np " + str(ncpus) + " python ../../main/stage3x_mpi_run.py 1000\n"
        print run_line
        fout.write(run_line)
    else:
        run_line = "mpirun -np " + str(ncpus) + " python ../../main/stage4x_mpi_run.py 5000\n"
        print run_line
        fout.write(run_line)

fout.close()
exit()
