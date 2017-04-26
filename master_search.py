#!/usr/bin/env python

"""
Project_Name: main, File_name: stage_1_search.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 14/04/15 , Time:01:05 PM

stage 1 in parallel
"""

import filters.pcs.pcsfilter as Pfilter
import filters.rdc.rdcfilter as Rfilter
import filters.noe.noefilter as Nfilter
import filters.noe.allNoes   as Noe
import filters.rmsd.RefRmsd as ref
import filters.sequence.sequence_similarity as Sfilter
import utility.io_util as io
import utility.smotif_util as sm
import utility.masterutil as mutil
import ranking.RDCStageRank as rank


def SmotifSearch(task):
    """
    Main ()
    :param index_array:
    :return:
    """

    index_array = task[0]
    stage = task[1]
    print index_array, stage
    s1_def, s2_def = mutil.getSSdef(index_array)
    smotif_def = sm.getSmotif(s1_def, s2_def)

    exp_data = io.readPickle("exp_data.pickle")
    exp_data_types = exp_data.keys()  # ['ss_seq', 'pcs_data', 'aa_seq', 'contacts']

    smotif_data = sm.readSmotifDatabase(smotif_def, exp_data['database_cutoff'])

    if not smotif_data:
        # If the smotif library doesn't exist, terminate further execution.
        return True

    dump_log = []
    # ************************************************************************************************
    # Main
    # The 'for' loop below iterates over all of the Smotifs and applies various filters
    # This is the place to add new filters as you desire. For starters, look at Sequence filter.
    # ************************************************************************************************


    for i in range(0, len(smotif_data)):

        # loop over for all of the entries in the smotif_db file

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

        if 'homologs' in exp_data_types:
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

        tlog, pcs_tensor_fits, rdc_tensor_fits, noe_fmeasure, = [], [], [], []
        ref_rmsd = 0.0
        tlog.append(['smotif', smotif_data[i]])
        tlog.append(['smotif_def', [s1_def, s2_def]])
        tlog.append(['cathcodes', [smotif_data[i][0]]])

        # ************************************************
        # Sequence filter
        # Aligns the smotif seq to target seq and calculates
        # sequence identity and the alignment score
        # ************************************************

        smotif_seq, seq_identity, blosum62_score = \
            Sfilter.SequenceSimilarity(s1_def, s2_def, smotif_data[i], exp_data)
        tlog.append(['seq_filter', smotif_seq, seq_identity, blosum62_score])

        # ************************************************
        # Unambiguous NOE score filter
        # uses experimental ambiguous noe data to filter Smotifs
        # scoring based on f-measure?
        # ************************************************

        if 'noe_data' in exp_data_types:
            # noe_fmeasure = Nfilter.s1NOEfit(s1_def, s2_def, smotif_data[i], exp_data)
            noe_fmeasure, no_of_noes = Noe.s1NOEfit(s1_def, s2_def, smotif_data[i], exp_data)
            tlog.append(['NOE_filter', noe_fmeasure, no_of_noes])
            tlog.append(['GlobalNoe_filter', noe_fmeasure, no_of_noes])

        # ************************************************
        # Residual dipolar coupling filter
        # uses experimental RDC data to filter Smotifs
        # scoring based on normalised chisqr.
        # ************************************************

        if 'rdc_data' in exp_data_types:
            rdc_tensor_fits = Rfilter.RDCAxRhFit(s1_def, s2_def, smotif_data[i], exp_data)
            if rdc_tensor_fits:
                tlog.append(['RDC_filter', rdc_tensor_fits])
            else:
                continue

        """
        if 'rdc_data' in exp_data_types:

            if noe_fmeasure and noe_fmeasure >= exp_data['noe_fmeasure'][0]:
                rdc_tensor_fits = Rfilter.RDCAxRhFit(s1_def, s2_def, smotif_data[i], exp_data)
                tlog.append(['RDC_filter', rdc_tensor_fits])

            elif seq_identity >=100:
                rdc_tensor_fits = Rfilter.RDCAxRhFit(s1_def, s2_def, smotif_data[i], exp_data)
                tlog.append(['RDC_filter', rdc_tensor_fits])
            else:
                continue
        """
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
        # if pcs_tensor_fits or rdc_tensor_fits:
        if pcs_tensor_fits:
            dump_log.append(tlog)

    # Save all of the hits in pickled arrays
    if dump_log:
        # testing for rdc_plus_pcs

        if 'rank_top_hits' in exp_data_types:
            dump_log = rank.rank_dump_log(dump_log, exp_data, stage=1)

        print "num of hits", len(dump_log)
        io.dumpPickle('0_' + str(index_array[0]) + "_" + str(index_array[1]) + ".pickle", dump_log)

    return True
