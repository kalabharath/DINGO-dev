#!/usr/bin/env python

"""
Project_Name: main, File_name: smotif_search.py
Aufthor: kalabharath, Email: kalabharath@gmail.com

"""
import filters.constraints.looplengthConstraint as llc
import filters.ilvaNOE.ilvanoepdf as noepdf
import filters.rdc.rdcfilter as Rfilter
import filters.rmsd.RefRmsd as ref
import filters.rmsd.qcp as qcp
import ranking.NoeStageRank as rank
import utility.io_util as io
import utility.alt_smotif_util as alt


def perform_alt_search(job, pair):
    # send_job = [tasks[t_job[0]], alt_sse_profile[t_job[1]], args.stage, task_index, lowest_noe_energy]

    dump_log = []
    task = job[0]
    ss_profile = job[1]
    old_noe_energy = job[-1]
    old_noe_energy = round(old_noe_energy, 3)
    stage = job[2]
    alt_smotif_log = task[9][1]
    direction = alt_smotif_log[-1]
    refine_pairs, computed_pairs = task[8][1], task[8][2]
    old_rdc_energy = task[6][3]
    old_rdc_energy = round(old_rdc_energy, 3)
    old_cath_codes = task[3][1]
    old_rmsd = task[7][1]
    smotif_coors, sse_ordered, rmsd = task[2][1], task[2][2], task[2][3]

    exp_data = io.getExpData()
    exp_data_types = exp_data.keys()  # ['ss_seq', 'pcs_data', 'aa_seq', 'contacts']

    # Check whether there are any noes for this pair
    print "Checking whether NOEs exist in this pair:", pair,
    if noepdf.noe_in_pair(sse_ordered, exp_data, pair):
        print "True"
        pass
    else:
        print "False"
        return False

    csmotif_data, sse_ordered, smotif_def = alt.getSmotifDB(sse_ordered, ss_profile, alt_smotif_log, pair,
                                                exp_data['database_cutoff'])

    if 'rmsd_cutoff' in exp_data_types:
        rmsd_cutoff = exp_data['rmsd_cutoff'][stage - 1]
    else:
        rmsd_cutoff = alt.altRMSDcutoff(smotif_def)

    if not csmotif_data:
        # If the smotif library doesn't exist.
        # Terminate further execution.
        return False

    # ************************************************************************************************
    # Main
    # The 'for' loop below iterates over all of the Smotifs and applies various filters
    # This is the place to add new filters as you desire. For starters, look at Sequence filter.
    # ************************************************************************************************

    for i in range(0, len(csmotif_data)):

        # ************************************************
        # Applying different filters for the Smotif assembly
        # ************************************************

        # Exclude natives if needed
        ref_rmsd, noe_probability = 0.0, 0.0
        tpdbid = csmotif_data[i][0][0]
        pdbid = tpdbid[0:4]

        if 'natives' in exp_data_types:
            natives = exp_data['natives']
            if pdbid in natives:
                continue
            # Stop further execution, but, iterate.
            else:
                pass

        if 'homologs' in exp_data_types:  # Smotif assembly only from the specified pdb files
            homologs = exp_data['homologs']
            if pdbid not in homologs:
                # Stop further execution, but, iterate.
                continue
            else:
                pass

        # ************************************************
        # RMSD filter using QCP method
        # quickly filters non-overlapping smotifs
        # ************************************************

        rmsd, transformed_coors = qcp.rmsdQCP4(pair, smotif_coors, alt_smotif_log, csmotif_data[i], direction,
                                               rmsd_cutoff)

        if rmsd <= rmsd_cutoff:

            # Loop constraint restricts the overlapping smotifs is not drifted far away.
            # loop_constraint = llc.loopConstraintAlt(transformed_coors, sse_ordered, direction)
            loop_constraint = True
            if loop_constraint:
                # Check whether the SSEs with in the assembled smotifs are clashing to one another
                no_clashes = qcp.kClahsesAltSmotif(transformed_coors, sse_ordered, pair, alt_smotif_log)
            else:
                no_clashes = False
        else:
            continue

        if no_clashes:

            # Prepare temporary arrays to log the data.
            tlog, total_percent, pcs_tensor_fits, rdc_tensor_fits = [], [], [], []
            tlog.append(['smotif', tpdbid])
            tlog.append(['smotif_def', sse_ordered])
            tlog.append(['qcp_rmsd', transformed_coors, sse_ordered, rmsd])
            tlog.append(['cathcodes', old_cath_codes])

            # ************************************************
            # Sequence filter
            # Aligns the smotif seq to target seq and calculates
            # sequence identity and the alignment score
            # ************************************************

            # concat current to previous seq
            # seq, seq_id = alt.getSeq(transformed_coors, sse_ordered, exp_data['aa_seq'])

            seq = 'SeqAnchor'
            seq_id = 30.0
            rmsd_cutoff = 1.0
            tlog.append(['seq_filter', seq, seq_id,exp_data['cluster_rmsd_cutoff']])

            # ************************************************
            # NOE score filter
            # uses experimental noe data to filter Smotifs
            # scoring based on log-likelihood?
            # ************************************************

            if 'ilva_noes' in exp_data_types:
                noe_probability, no_of_noes, noe_energy, noe_data, new_cluster_protons, new_cluster_sidechains = noepdf.refineILVA(
                    transformed_coors, sse_ordered, exp_data, old_noe_energy)

                if noe_probability >= exp_data['expected_noe_prob'][stage - 1]:
                    noe_energy = round(noe_energy, 3)
                    tlog.append(['NOE_filter', noe_probability, no_of_noes, noe_energy, noe_data, new_cluster_protons,
                                 new_cluster_sidechains])
                else:
                    continue

            # ************************************************
            # Residual dipolar coupling filter
            # uses experimental RDC data to filter Smotifs
            # scoring based on normalised chisqr
            # ************************************************

            if 'rdc_data' in exp_data_types:
                rdc_tensor_fits, log_likelihood, rdc_energy = Rfilter.RDCAxRhFit2(transformed_coors, sse_ordered,
                                                                                  exp_data, stage)

                if rdc_energy == 999.99:
                    continue
                else:
                    tlog.append(['RDC_filter', rdc_tensor_fits, log_likelihood, rdc_energy])

            # ************************************************
            # Calc RMSD of the reference structure.
            # Used to identify the lowest possible RMSD
            # structure for the target, from the Smotif library.
            # ************************************************

            if 'reference_ca' in exp_data_types:
                ref_rmsd = ref.calcRefRMSD2(exp_data['reference_ca'], sse_ordered, transformed_coors)
                tlog.append(['Ref_RMSD', ref_rmsd, seq_id])
                log_refine_pair = [pair, tpdbid]
                try:
                    log_refine_smotif = (task[8][3])[:]
                except:
                    log_refine_smotif = []

                log_refine_smotif.append(log_refine_pair)
                tlog.append(['Refine_smotifs', refine_pairs, computed_pairs, log_refine_smotif])

            if (noe_energy < old_noe_energy) or (rdc_energy < old_rdc_energy):
                print "rmsd:", rmsd, pair
                print "NOE energy", old_noe_energy, noe_energy, noe_probability
                print "RDC energy", old_rdc_energy, rdc_energy
                print "Ref_rmsd", old_rmsd, ref_rmsd
                dump_log.append(tlog)
            elif (noe_energy == old_noe_energy) and (rdc_energy < old_rdc_energy):
                print "rmsd:", rmsd, pair
                print "NOE energy", old_noe_energy, noe_energy, noe_probability
                print "RDC energy", old_rdc_energy, rdc_energy
                print "Ref_rmsd", old_rmsd, ref_rmsd
                dump_log.append(tlog)
            elif (rdc_energy == old_rdc_energy) and (noe_energy < old_noe_energy):
                print "rmsd:", rmsd, pair
                print "NOE energy", old_noe_energy, noe_energy, noe_probability
                print "RDC energy", old_rdc_energy, rdc_energy
                print "Ref_rmsd", old_rmsd, ref_rmsd
                dump_log.append(tlog)
            else:
                continue

    # Dumping hits as a pickle array.
    if len(dump_log) > 0:
        if 'rank_top_hits' in exp_data_types:
            rank_top_hits = exp_data['rank_top_hits']
            num_hits = rank_top_hits[stage - 1]
            dump_log = rank.rank_assembly_with_clustering(dump_log, num_hits)

            print "Reducing the amount of data to:", rank_top_hits[stage - 1], len(dump_log)
        print "num of hits", len(dump_log),

        return dump_log
    else:
        return False


def altSmotifSearch(job):

    # send_job = [tasks[t_job[0]], alt_sse_profile[t_job[1]], args.stage, task_index, lowest_noe_energy]
    all_log = []
    task = (job[0])[:]
    refine_pair = task[8][1]
    index_array = task[3]
    print index_array
    for pair in refine_pair:
        tdump_log = perform_alt_search(job, pair)
        if tdump_log:
            for t in tdump_log:
                all_log.append(t)

    # Dump data to the disk
    if all_log:
        io.dumpGzipPickle("rtx_" + str(index_array[0]) + "_" + str(index_array[1]) + ".gzip", all_log)
        return True
    else:
        return False
