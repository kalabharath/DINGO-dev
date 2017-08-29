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
import utility.PCSmap     as PCSmap
import utility.referenceUtil as ref

# Parse the input data text file
data = io.readInputDataFiles('input_data.txt')
datatypes = data.keys()

data_dict = {}
print datatypes

if 'fasta_file' in datatypes:
    # Read in the amino acid sequence and the Secondary structure assignment
    handle, aa_seq = io.readFasta(data['fasta_file'])
    ss_seq = ru.matchSeq2SS(aa_seq, data['ss_file'])
    # Generate fuzzy +/-2 SSE combinations
    ss_def, ss_combi = ss.genSSCombinations(ss_seq)
    io.dumpPickle("ss_profiles.pickle", ss_combi)
    data_dict['ss_seq'] = ss_seq
    data_dict['aa_seq'] = aa_seq
    print ss_def
    print ss_seq
else:
    pass

if 'native_pdbs' in datatypes:
    # Read the native pdbs that you can exclude from the smotif search
    native_pdbs = data['native_pdbs']
    native_pdbs = native_pdbs.lower()
    native_pdbs = native_pdbs.split()
    data_dict['natives'] = native_pdbs
    print native_pdbs
else:
    pass

if 'homolog_pdbs' in datatypes:
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
    data_dict['homologs'] = thomologs
else:
    pass

if 'predicted_axial' in datatypes:
    pred_axial = data['predicted_axial']
    pred_axial = pred_axial.split()
    pred_axial = [float(i) for i in pred_axial]
    data_dict['pred_axial'] = pred_axial
else:
    pass

if 'TinK' in datatypes:
    TinK = data['TinK']
    TinK = TinK.split()
    TinK = [float(i) for i in TinK]
    data_dict['TinK'] = TinK
    print TinK
else:
    pass

if 'B0inT' in datatypes:
    B0inT = data['B0inT']
    B0inT = B0inT.split()
    B0inT = [float(i) for i in B0inT]
    data_dict['B0inT'] = B0inT
    print B0inT
else:
    pass

if 'exp_error' in datatypes:
    exp_error = data['exp_error']
    exp_error = exp_error.split()
    exp_error = [float(i) for i in exp_error]
    data_dict['exp_error'] = exp_error
else:
    pass

if 'abs_exp_error' in datatypes:
    abs_exp_error = data['abs_exp_error']
    abs_exp_error = abs_exp_error.split()
    abs_exp_error = [float(i) for i in abs_exp_error]
    data_dict['abs_exp_error'] = abs_exp_error
else:
    pass

if 'rank_top_hits' in datatypes:
    rank_top_hits = data['rank_top_hits']
    rank_top_hits = rank_top_hits.split()
    rank_top_hits = [float(i) for i in rank_top_hits]
    data_dict['rank_top_hits'] = rank_top_hits
else:
    pass

if 'noe_fmeasure' in datatypes:
    noe_fmeasure = data['noe_fmeasure']
    noe_fmeasure = noe_fmeasure.split()
    noe_fmeasure = [float(i) for i in noe_fmeasure]
    print "noe_fmeasure", noe_fmeasure
    data_dict['noe_fmeasure'] = noe_fmeasure
else:
    pass

if 'rmsd_cutoff' in datatypes:
    rmsd_cutoff = data['rmsd_cutoff']
    rmsd_cutoff = rmsd_cutoff.split()
    rmsd_cutoff = [float(i) for i in rmsd_cutoff]
    data_dict['rmsd_cutoff'] = rmsd_cutoff
else:
    pass

if 'reference_pdb' in datatypes:
    reference_ca = ref.getRefCoors(data['reference_pdb'])
    data_dict['reference_ca'] = reference_ca
    print rmsd_cutoff
else:
    pass

if 'clash_distance' in datatypes:
    clash_distance = float(data['clash_distance'])
    data_dict['clash_distance'] = clash_distance
    print 'clash_distance: ', clash_distance

if 'rdc_input_files' in datatypes:
    rdc_data = ru.getRdcData(data['rdc_input_files'], ss_seq)
    data_dict['rdc_data'] = rdc_data
else:
    pass

if 'noe_bmrb_mr' in datatypes:
    if 'bmrb_sequence_correction' in datatypes:
        seq_correction = int(data['bmrb_sequence_correction'])
        noe_data, total_noe_count = nu.parseBMRBblockMR(data['noe_bmrb_mr'], data_dict['aa_seq'],seq_correction)
        data_dict['noe_data'] = [noe_data, total_noe_count]
else:
    pass

if 'noe_input_files' in datatypes:
    # noe_data, total_noe_count = nu.getNOEData(data['noe_input_files'], ss_seq)
    noe_data, total_noe_count = nu.parseNOEData(data['noe_input_files'])
    data_dict['noe_data'] = [noe_data, total_noe_count]
else:
    pass

if 'pcs_broker' in datatypes:
    print data['pcs_broker']
    pcs_broker = data['pcs_broker']
    pcsdata = io.getPcsTagInfo(ss_seq, pcs_broker)
    data_dict['pcs_data'] = pcsdata
else:
    pass

if 'axrh_cutoff' in datatypes:
    axrh_cutoff = data['axrh_cutoff']
    axrh_cutoff = axrh_cutoff.split()
    axrh_cutoff = [float(i) for i in axrh_cutoff]
    data_dict['axrh_cutoff'] = axrh_cutoff

if 'chisqr_cutoff' in datatypes:
    chisqr_cutoff = data['chisqr_cutoff']
    chisqr_cutoff = chisqr_cutoff.split()
    chisqr_cutoff = [float(i) for i in chisqr_cutoff]
    data_dict['chisqr_cutoff'] = chisqr_cutoff

if ('rdc_input_files' in datatypes) and ('pcs_broker' in datatypes):
    map_route = PCSmap.getRoute(ss_seq, pcsdata)
    print map_route
    io.dumpPickle("pcs_route.pickle", map_route)
elif ('pcs_broker' in datatypes) and not ('rdc_input_files' in datatypes):
    map_route = PCSmap.getRoute(ss_seq, pcsdata)
    print map_route
    io.dumpPickle("pcs_route.pickle", map_route)
else:
    map_route = ru.getRDCMapRoute(ss_combi, rdc_data)
    print map_route
    io.dumpPickle("rdc_route.pickle", map_route)

if 'database_cutoff' in datatypes:
    database_cutoff = data['database_cutoff']
    data_dict['database_cutoff'] = database_cutoff
else:
    pass

if 'metal_spheres' in datatypes:
    metals_def = data['metal_spheres']
    metals_def = metals_def.split()
    metals_def = [int(i) for i in metals_def]
    print metals_def
    npts = PCSmap.countPointsOnSpheres([metals_def[0], metals_def[1]], metals_def[2])
    metals_def.append(npts)
    print metals_def
    data_dict['metal_spheres'] = metals_def
else:
    # default metal definitions
    metals_def = [1, 45, 200, 198000]
    data_dict['metal_spheres'] = metals_def
    pass

io.dumpPickle("exp_data.pickle", data_dict)

fout = open("run.sh", 'w')
fout.write("#!/bin/bash\n")
ncpus = multiprocessing.cpu_count()

for i in range(0, len(map_route)):
    smotif = map_route[i]
    print smotif
    if i == 0:
        run_line = "mpirun -np " + str(ncpus) + " python ../../main/Dingo_mpi.py --stage 1" + "  --numhits 0 \n"
        print run_line
        fout.write(run_line)
    elif i == 1:
        run_line = "mpirun -np " + str(ncpus) + " python ../../main/Dingo_mpi.py --stage 2" + "  --numhits 1000 \n"
        print run_line
        fout.write(run_line)
    elif i != 1 and i <= len(map_route) - 3:
        run_line = "mpirun -np " + str(ncpus) + " python ../../main/Dingo_mpi.py --stage 3" + "  --numhits 1000 \n"
        print run_line
        fout.write(run_line)
    else:
        run_line = "mpirun -np " + str(ncpus) + " python ../../main/Dingo_mpi.py --stage 4" + "  --numhits 1000 \n"
        print run_line
        fout.write(run_line)
fout.close()
exit()
