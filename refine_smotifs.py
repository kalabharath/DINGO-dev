import utility.io_util as io
import filters.rmsd.refine_rmsd as qcp
import filters.ilvaNOE.ilvanoepdf as noepdf
import filters.rdc.rdcfilter as Rfilter
import filters.rmsd.RefRmsd as ref
import utility.smotif_util as sm


def getRefinementIndices(sse_array):
    import itertools
    indices = list(itertools.combinations(range(len(sse_array)), 2))
    refine_pairs = []
    for pair in indices:
        if abs(pair[0] - pair[1]) == 1:
            pass
        else:
            t_array = sm.array2string([sse_array[pair[0]], sse_array[pair[1]]])
            # print t_array
            if t_array in refine_pairs:
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

def getSeq(coor_array, sse_ordered, aa_seq):

    one_letter = {'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'GLU': 'E', 'GLN': 'Q',
                  'ASP': 'D', 'ASN': 'N', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F', 'TYR': 'Y',
                  'ARG': 'R', 'LYS': 'K', 'SER': 'S', 'THR': 'T', 'MET': 'M', 'ALA': 'A',
                  'GLY': 'G', 'PRO': 'P', 'CYS': 'C', 'ASX': 'D', 'GLX': 'G', 'UNK': 'A'}
    concat_seq = ''
    for frag in coor_array:
        atom_num = 1
        for i in range(atom_num, len(frag[0]), 5):
            res = (frag[5][i])
            concat_seq = concat_seq+one_letter[res]

    native_sse_seq = ''
    for sse in sse_ordered:
        sse_seq = aa_seq[sse[4] - 1: sse[5]]
        native_sse_seq = native_sse_seq + sse_seq
    k = 0.0
    if len(concat_seq) != len(native_sse_seq):
        print "Something is wrong with extracting sequence information"
    for i in range(0, len(concat_seq)):

        if native_sse_seq[i] == concat_seq[i]:
            k += 1
    seq_id = (k / float(len(concat_seq))) * 100
    return concat_seq, seq_id

def performRefinement(task, stage, pair):

    #task = work[0]
    #stage = work[1]
    #task_index = work[2]
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
    8: Refine_Smotifs

    """

    smotif_coors, sse_ordered, rmsd = task[2][1], task[2][2], task[2][3]
    refine_pairs, computed_pairs = task[8][1], task[8][2]
    old_noe_energy = task[5][3]
    old_rdc_energy = task[6][3]
    old_cath_codes, parent_smotifs = task[3][1], task[3][2]
    old_rmsd = task[7][1]
    old_noe_energy = round(old_noe_energy, 3)

    tdump_log = []

    if old_noe_energy <= 0.005:
        print "NOE energy is Zero there is no need to do any refinement, exiting task:"
        return tdump_log
    else:
        print "Energy is nonzero proceeding with refinement: ", old_noe_energy


    db_entries = getfromDB(pair, sse_ordered, exp_data['database_cutoff'])

    for smotif in db_entries:
        tlog = []
        rmsd_cutoff = exp_data['rmsd_cutoff'][stage - 1]
        transformed_coors, rmsd = qcp.refineRMSD(smotif_coors, pair, smotif, rmsd_cutoff)
        if rmsd <= rmsd_cutoff:
            if not qcp.kClahsesRefined(transformed_coors, sse_ordered, pair):
                continue
        else:
            continue

        seq, seq_id = getSeq(transformed_coors, sse_ordered, exp_data['aa_seq'])
        print seq,seq_id
        die

        tlog.append(['smotif', smotif])
        tlog.append(['smotif_def', sse_ordered])
        tlog.append(['qcp_rmsd', transformed_coors, sse_ordered, rmsd])
        tlog.append(['cathcodes', old_cath_codes, parent_smotifs])
        tlog.append(['seq_filter',seq, seq_id ])

        # Recalculate NOE energy
        if 'ilva_noes' in exp_data_types:
            noe_probability, no_of_noes, noe_energy, noe_data, new_cluster_protons, \
            new_cluster_sidechains = noepdf.refineILVA(transformed_coors, sse_ordered, exp_data, stage)

            if noe_probability >= exp_data['expected_noe_prob'][stage - 1]:
                if noe_energy <= old_noe_energy:
                    tlog.append(
                        ['NOE_filter', noe_probability, no_of_noes, noe_energy, noe_data, new_cluster_protons,
                         new_cluster_sidechains])
                else:
                    continue
            else:
                continue

        if 'rdc_data' in exp_data_types:
            rdc_tensor_fits, log_likelihood, rdc_energy = Rfilter.RDCAxRhFit2(transformed_coors, sse_ordered,
                                                                              exp_data, stage)
            tlog.append(['RDC_filter', rdc_tensor_fits, log_likelihood, rdc_energy])

        if 'reference_ca' in exp_data_types:
            ref_rmsd = ref.calcRefRMSD2(exp_data['reference_ca'], sse_ordered, transformed_coors)
            tlog.append(['Ref_RMSD', ref_rmsd, seq_id])
            tlog.append(['Refine_Smotifs', refine_pairs, computed_pairs])
            print "rmsd:", rmsd, pair
            print "NOE energy", old_noe_energy, noe_energy, noe_probability
            print "RDC energy", old_rdc_energy, rdc_energy
            print "Ref_rmsd", old_rmsd, ref_rmsd

        tdump_log.append(tlog)
        if len(tdump_log) >= 5:
            break

        # copy the new coordinates for the next sse pair
    if tdump_log:

        return tdump_log
    else:
        return False


def SmotifRefinement(work):
    task = work[0]
    stage = work[1]
    task_index = work[2]


    #smotif_coors, sse_ordered, rmsd = task[2][1], task[2][2], task[2][3]

    refine_pairs = task[8][1]
    print refine_pairs

    #refine_pairs = getRefinementIndices(sse_ordered)
    old_noe_energy = task[5][3]
    old_noe_energy = round(old_noe_energy, 3)


    dump_log = []


    if old_noe_energy <= 0.005:
        print "NOE energy is Zero there is no need to do any refinement, exiting task:", task_index

        dump_log.append(task)
        return dump_log
    else:
        dump_log.append(task)
        print "Energy is nonzero proceeding with refinement: ", old_noe_energy

    print "All pairs: ", refine_pairs

    refine_pairs =  [(0, 2), (0, 3)]

    for pair in refine_pairs:
        t_log = []
        for entry in dump_log:
            print "testing start", len(dump_log), pair
            t_entry = performRefinement(entry, stage, pair)
            if t_entry:
                for t in t_entry:
                    t_log.append(t)
        for t in t_log:
            dump_log.append(t)

    return dump_log
