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


def S1SmotifSearch(task):
    """
    Main ()
    :param task:
    :return:
    """

    index_array = task[0]
    stage = task[1]
    s1_def, s2_def = mutil.getSSdef(index_array)
    smotif_def = sm.getSmotif(s1_def, s2_def)

    exp_data = io.readPickle("exp_data.pickle")
    exp_data_types = exp_data.keys()  # ['ss_seq', 'pcs_data', 'aa_seq', 'contacts']

    smotif_data = sm.readSmotifDatabase(smotif_def, exp_data['database_cutoff'])

    if not smotif_data:
        # If the smotif library doesn't exist, terminate further execution.
        return False

    dump_log = []

    # ************************************************************************************************
    # Main
    # The 'for' loop below iterates over all of the Smotifs and applies various filters
    # This is the place to add new filters as you desire. For starters, look at Sequence filter.
    # ************************************************************************************************


    for i in range(0, len(smotif_data)):

        # ************************************************
        # Excluding the natives
        # ************************************************

        natives = exp_data['natives']
        tpdbid = smotif_data[i][0][0]
        pdbid = tpdbid[0:4]

        if 'natives' in exp_data_types:
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
        # Applying different filters to Smotifs
        # Prepare temp log array to save data at the end
        # ************************************************
        # initialize variables
        tlog, pcs_tensor_fits, rdc_tensor_fits, = [], [], []
        ref_rmsd, noe_probability = 0.0, 0.0

        tlog.append(['smotif', smotif_data[i]])
        tlog.append(['smotif_def', [s1_def, s2_def]])
        parent_smotifs = sm.array2string([smotif_data[i][0]])
        tlog.append(['qcp_rmsd'])
        tlog.append(['cathcodes', [smotif_data[i][0]], parent_smotifs])
        # ************************************************
        # Sequence filter
        # Aligns the smotif seq to target seq and calculates
        # sequence identity and the alignment score
        # ************************************************
        smotif_seq, seq_identity = Sfilter.getS1SeqIdentity(s1_def, s2_def, smotif_data[i], exp_data)
        tlog.append(['seq_filter', smotif_seq, seq_identity])

        # ************************************************
        # Unambiguous NOE score filter
        # uses experimental ambiguous noe data to filter Smotifs
        # scoring based on f-measure?
        # ************************************************

        if 'ilva_noes' in exp_data_types:
            noe_probability, no_of_noes, noe_energy, noe_data, cluster_protons, cluster_sidechains = noepdf.s1ILVApdf(s1_def, s2_def, smotif_data[i], exp_data, stage)

            if noe_probability >= exp_data['expected_noe_prob'][stage - 1]:
                tlog.append(['NOE_filter', noe_probability, no_of_noes, noe_energy, noe_data, cluster_protons, cluster_sidechains])
            else:
                continue

        # ************************************************
        # Residual dipolar coupling filter
        # uses experimental RDC data to filter Smotifs
        # scoring based on normalised chisqr.
        # ************************************************

        if 'rdc_data' in exp_data_types:
            rdc_tensor_fits, log_likelihood, rdc_energy = Rfilter.RDCAxRhFit(s1_def, s2_def, smotif_data[i], exp_data)
            if rdc_tensor_fits:
                tlog.append(['RDC_filter', rdc_tensor_fits, log_likelihood, rdc_energy])
            else:
                continue

        # ************************************************
        # Pseudocontact Shift filter
        # uses experimental PCS data to filter Smotifs
        # scoring based on normalised chisqr
        # ************************************************

        if 'pcs_data' in exp_data_types:
            pcs_tensor_fits = Pfilter.PCSAxRhFit(s1_def, s2_def, smotif_data[i], exp_data)
            tlog.append(['PCS_filter', pcs_tensor_fits])

        # ************************************************
        # Calc RMSD of the reference structure.
        # Used to identify the lowest possible RMSD
        # structure for the target, from the Smotif library.
        # ************************************************

        if 'reference_ca' in exp_data_types:
            ref_rmsd = ref.calcRefRMSD(exp_data['reference_ca'], s1_def, s2_def, smotif_data[i], rmsd_cutoff=100.0)
            tlog.append(['Ref_RMSD', ref_rmsd, seq_identity])

        # Dump the data to the disk
        if pcs_tensor_fits or noe_probability:
            dump_log.append(tlog)

    # Save all of the hits in pickled arrays
    if dump_log:
        # testing for rdc_plus_pcs
        if 'rank_top_hits' in exp_data_types:
            #dump_log = rank.rank_dump_log(dump_log, exp_data, stage=1)
            rank_top_hits = exp_data['rank_top_hits']
            num_hits = rank_top_hits[stage - 1]
            dump_log = rank.rank_assembly(dump_log, num_hits)
            print "Reducing the amount of data to:", rank_top_hits[stage - 1], len(dump_log)
        print "num of hits", len(dump_log)
        # io.dumpPickle('0_' + str(index_array[0]) + "_" + str(index_array[1]) + ".pickle", dump_log)
        io.dumpGzipPickle('0_' + str(index_array[0]) + "_" + str(index_array[1]) + ".gzip", dump_log)
        return dump_log
    else:
        return False


def sXSmotifSearch(task):
    """
     Main()
    :param task:
    :return:
    """
    index_array = task[0]
    stage = task[1]

    exp_data = io.readPickle("exp_data.pickle")
    exp_data_types = exp_data.keys()  # ['ss_seq', 'pcs_data', 'aa_seq', 'contacts']
    psmotif, preSSE, dump_log = [], [], []

    if stage == 2:
        psmotif = uts2.getPreviousSmotif(index_array[0])
        current_ss, direction = uts2.getSS2(index_array[1])
        csmotif_data, smotif_def = mutil.getfromDB(psmotif, current_ss, direction, exp_data['database_cutoff'], stage)
        sse_ordered, refine_pairs, computed_pairs = mutil.orderSSE(psmotif, current_ss, direction, stage)
        sorted_noe_data, cluster_protons, cluster_sidechains = mutil.fetchNOEdata(psmotif)
    else:

        preSSE = uts2.getPreviousSmotif(index_array[0])
        current_ss, direction = uts2.getSS2(index_array[1])
        csmotif_data, smotif_def = mutil.getfromDB(preSSE, current_ss, direction, exp_data['database_cutoff'], stage)
        sse_ordered, refine_pairs, computed_pairs = mutil.orderSSE(preSSE, current_ss, direction, stage)
        sorted_noe_data, cluster_protons, cluster_sidechains = mutil.fetchNOEdata(preSSE)

    print current_ss, direction
    if not csmotif_data:
        # If the smotif library doesn't exist.
        # Terminate further execution.
        return True

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

        if stage == 2:
            rmsd, transformed_coos = qcp.rmsdQCP(psmotif[0], csmotif_data[i], direction, rmsd_cutoff)
        else:
            rmsd, transformed_coos = qcp.rmsdQCP3(preSSE, csmotif_data[i], direction, rmsd_cutoff)

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
            tlog.append(['smotif', csmotif_data[i]])
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
                noe_probability, no_of_noes, noe_energy, noe_data, new_cluster_protons, new_cluster_sidechains = noepdf.sX2ILVApdf(transformed_coos,
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
                rdc_tensor_fits, log_likelihood, rdc_energy = Rfilter.RDCAxRhFit2(transformed_coos, sse_ordered, exp_data, stage)
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
                tlog.append(['Refine_Smotifs', refine_pairs, computed_pairs])
                print tlog [-1]

            if pcs_tensor_fits or noe_probability:
                # dump data to the disk
                dump_log.append(tlog)

    # Dumping hits as a pickle array.
    if len(dump_log) > 0:
        if 'rank_top_hits' in exp_data_types:
            #dump_log = rank.rank_dump_log(dump_log, exp_data, stage)
            rank_top_hits = exp_data['rank_top_hits']
            num_hits = rank_top_hits[stage - 1]
            dump_log = rank.rank_assembly(dump_log, num_hits)
            print "Reducing the amount of data to:", rank_top_hits[stage - 1], len(dump_log)
        print "num of hits", len(dump_log),
        # io.dumpPickle("tx_" + str(index_array[0]) + "_" + str(index_array[1]) + ".pickle", dump_log)
        #io.dumpGzipPickle("tx_" + str(index_array[0]) + "_" + str(index_array[1]) + ".gzip", dump_log)
        return dump_log
    else:
        return False
