#!/usr/bin/env python

"""
Project_Name: main, File_name: smotif_search.py
Aufthor: kalabharath, Email: kalabharath@gmail.com

"""
import filters.constraints.looplengthConstraint as llc
import filters.ilvaNOE.ilvanoepdf as noepdf
import filters.pcs.pcsfilter as Pfilter
import filters.rdc.rdcfilter as Rfilter
import filters.rmsd.RefRmsd as ref
import filters.rmsd.qcp as qcp
import filters.sequence.sequence_similarity as Sfilter
import ranking.NoeStageRank as rank
import utility.io_util as io
import utility.masterutil as mutil
import utility.smotif_util as sm
import utility.stage2_util as uts2
import utility.alt_smotif_util as alt


def perform_alt_search(job, pair):
    # send_job = [tasks[t_job[0]], alt_sse_profile[t_job[1]], args.stage, task_index, lowest_noe_energy]

    dump_log = []
    task = job[0]
    ss_profile = job[1]
    old_noe_energy = job[-1]
    stage = job[2]
    alt_smotif_log = task[9][1]
    direction = alt_smotif_log[-1]

    exp_data = io.getExpData()
    exp_data_types = exp_data.keys()  # ['ss_seq', 'pcs_data', 'aa_seq', 'contacts']

    smotif_coors, sse_ordered, rmsd = task[2][1], task[2][2], task[2][3]
    print len(smotif_coors), sse_ordered, rmsd
    psmotif, preSSE, = [], []

    # Check whether there are any noes for this pair
    if noepdf.noe_in_pair(sse_ordered, exp_data, pair):
        pass
    else:
        return False

    csmotif_data = alt.getSmotifDB(sse_ordered, ss_profile, alt_smotif_log, pair, exp_data['database_cutoff'])
    print "Successfully parsed db",len(csmotif_data)

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
        no_clashes = False

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

        rmsd_cutoff = exp_data['rmsd_cutoff'][stage - 1]

        # preSSEs shoud be modified
        trunk_sse_coors = alt.delete_last_sse(smotif_coors, alt_smotif_log)
        print "len trunk_sse_coors", len(trunk_sse_coors)
        print csmotif_data[i][0]
        rmsd, transformed_coos = qcp.rmsdQCP4(pair, trunk_sse_coors, csmotif_data[i], direction, rmsd_cutoff)
        print "New RMSD", rmsd
        
        if rmsd <= rmsd_cutoff:
            # Loop constraint restricts the overlapping smotifs is not drifted far away.
            loop_constraint = llc.loopConstraint(transformed_coos, sse_ordered, direction, smotif_def)
            if loop_constraint:
                # Check whether the SSEs with in the assembled smotifs are clashing to one another
                no_clashes = qcp.kClashes(transformed_coos, sse_ordered, current_ss)
            else:
                no_clashes = False
        else:
            continue

        if no_clashes:

            # Prepare temporary arrays to log the data.
            tlog, total_percent, pcs_tensor_fits, rdc_tensor_fits = [], [], [], []
            tlog.append(['smotif', tpdbid])
            tlog.append(['smotif_def', sse_ordered])
            tlog.append(['qcp_rmsd', transformed_coos, sse_ordered, rmsd])

            if stage == 2:
                cathcodes, parent_smotifs = sm.orderCATH(psmotif, csmotif_data[i][0], direction)
            else:
                cathcodes, parent_smotifs = sm.orderCATH(preSSE, csmotif_data[i][0], direction)
            tlog.append(['cathcodes', cathcodes, parent_smotifs])

            # ************************************************
            # Sequence filter
            # Aligns the smotif seq to target seq and calculates
            # sequence identity and the alignment score
            # ************************************************

            # concat current to previous seq
            if stage == 2:
                seq_identity, concat_seq = Sfilter.getSXSeqIdentity(current_ss, csmotif_data[i], direction, exp_data,
                                                                    psmotif, sse_ordered)
            else:
                seq_identity, concat_seq = Sfilter.getSXSeqIdentity(current_ss, csmotif_data[i], direction, exp_data,
                                                                    preSSE, sse_ordered)

            tlog.append(['seq_filter', concat_seq, seq_identity])

            # ************************************************
            # NOE score filter
            # uses experimental noe data to filter Smotifs
            # scoring based on log-likelihood?
            # ************************************************

            if 'ilva_noes' in exp_data_types:
                noe_probability, no_of_noes, noe_energy, noe_data, new_cluster_protons, new_cluster_sidechains = noepdf.sX2ILVApdf(
                    transformed_coos,
                    sse_ordered, current_ss,
                    sorted_noe_data,
                    cluster_protons, cluster_sidechains, exp_data, stage)

                if noe_probability >= exp_data['expected_noe_prob'][stage - 1]:
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
                rdc_tensor_fits, log_likelihood, rdc_energy = Rfilter.RDCAxRhFit2(transformed_coos, sse_ordered,
                                                                                  exp_data, stage)
                if rdc_tensor_fits:
                    tlog.append(['RDC_filter', rdc_tensor_fits, log_likelihood, rdc_energy])
                else:
                    # Do not execute any further
                    continue

            # ************************************************
            # Pseudocontact Shift filter
            # uses experimental PCS data to filter Smotifs
            # scoring based on normalised chisqr
            # ************************************************

            if 'pcs_data' in exp_data_types:
                pcs_tensor_fits = Pfilter.PCSAxRhFit2(transformed_coos, sse_ordered, exp_data, stage)
                tlog.append(['PCS_filter', pcs_tensor_fits])

            # ************************************************
            # Calc RMSD of the reference structure.
            # Used to identify the lowest possible RMSD
            # structure for the target, from the Smotif library.
            # ************************************************

            if 'reference_ca' in exp_data_types:
                ref_rmsd = ref.calcRefRMSD2(exp_data['reference_ca'], sse_ordered, transformed_coos)
                tlog.append(['Ref_RMSD', ref_rmsd, seq_identity])
                tlog.append(['Refine_Smotifs', refine_pairs, computed_pairs, log_refine_smotif])

            if pcs_tensor_fits or noe_probability:
                # dump data to the disk
                dump_log.append(tlog)

    # Dumping hits as a pickle array.
    if len(dump_log) > 0:
        if 'rank_top_hits' in exp_data_types:
            rank_top_hits = exp_data['rank_top_hits']
            num_hits = rank_top_hits[stage - 1]
            dump_log = rank.rank_assembly_with_clustering(dump_log, exp_data['aa_seq'], num_hits)

            print "Reducing the amount of data to:", rank_top_hits[stage - 1], len(dump_log)
        print "num of hits", len(dump_log),

        return dump_log
    else:
        return False


def altSmotifSearch(job):
    # send_job = [tasks[t_job[0]], alt_sse_profile[t_job[1]], args.stage, task_index, lowest_noe_energy]

    tdump_log =[]
    task = job[0]
    smotif_refinement = task[8]
    refine_pair = smotif_refinement[1]
    print "refine_pair", refine_pair

    for pair in refine_pair:
        tdump_log = perform_alt_search(job, pair)

    tdump_log.append(task)