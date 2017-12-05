import utility.io_util as io
import filters.rmsd.refine_rmsd as qcp
import filters.ilvaNOE.ilvanoepdf as noepdf
import filters.rdc.rdcfilter as Rfilter


def getRefinementIndices(sse_array):
    import itertools
    indices = list(itertools.combinations(range(len(sse_array)), 2))
    refine_pairs = []
    for pair in indices:
        if abs(pair[0] - pair[1]) == 1:
            pass
        else:
            refine_pairs.append(pair)
    return refine_pairs


def getfromDB(pair, sse_ordered, database_cutoff):
    from utility.smotif_util import getSmotif, readSmotifDatabase
    s1 = sse_ordered[pair[0]]
    s2 = sse_ordered[pair[1]]
    smotif = getSmotif(s1, s2)
    return readSmotifDatabase(smotif, database_cutoff)


def SmotifRefinement(work):
    task = work[0]
    stage = work[1]
    task_index = work[2]
    exp_data = io.readPickle("exp_data.pickle")
    exp_data_types = exp_data.keys()  # ['ss_seq', 'pcs_data', 'aa_seq', 'contacts']

    """              
    0: smotif
    1: smotif_def
    2: qcp_rmsd
    3: cathcodes
    4: seq_filter
    5: NOE_filter
    6: RDC_filter
    7: Ref_RMSD   
    
    """

    smotif_coors, sse_ordered, rmsd = task[2][1], task[2][2], task[2][3]
    refine_pairs = getRefinementIndices(sse_ordered)
    old_noe_energy = task[5][3]
    old_rdc_energy = task[6]

    dump_log = []
    dump_log.append(task)

    for pair in refine_pairs:

        tlog = []

        db_entries = getfromDB(pair, sse_ordered, exp_data['database_cutoff'])
        for smotif in db_entries:
            rmsd_cutoff = exp_data['rmsd_cutoff'][stage - 1]
            transformed_coors, rmsd = qcp.refineRMSD(smotif_coors, pair, smotif, rmsd_cutoff)
            if rmsd <= rmsd_cutoff:

                if qcp.kClahsesRefined(transformed_coors, sse_ordered, pair):
                    pass
                else:
                    continue

                tlog.append(['smotif', smotif[0]])
                tlog.append(['smotif_def', sse_ordered])
                tlog.append(['qcp_rmsd', transformed_coors, sse_ordered, rmsd])
                tlog.append(['cathcodes','figure this out'])
                tlog.append(['seq_filter', "figure this out"])
                # check for clashes

                    # Recalculate NOE energy
                if 'ilva_noes' in exp_data_types:
                    noe_probability, no_of_noes, noe_energy, noe_data, new_cluster_protons,\
                    new_cluster_sidechains = noepdf.refineILVA(transformed_coors, sse_ordered, exp_data, stage)

                    if noe_probability >= exp_data['expected_noe_prob'][stage - 1]:
                        print "rmsd:", rmsd
                        print "NOE energy", old_noe_energy, noe_energy, noe_probability
                        if noe_energy <= old_noe_energy
                            tlog.append(['NOE_filter', noe_probability, no_of_noes, noe_energy, noe_data, new_cluster_protons, new_cluster_sidechains])
                        else:
                            continue

                if 'rdc_data' in exp_data_types:
                    rdc_tensor_fits, log_likelihood = Rfilter.RDCAxRhFit2(transformed_coors, sse_ordered,exp_data, stage)
                    if rdc_tensor_fits:
                        dump_log.append(tlog)


    io.dumpGzipPickle("tx_" + str(task_index) + ".gzip", dump_log)
    return False
