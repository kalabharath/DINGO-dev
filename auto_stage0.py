#!/usr/bin/env python

"""
Project_Name: main, File_name: stage_0.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 31/03/15 , Time:10:56 AM

Prepare all the relavant files for stage1 & 2

"""
import sys
sys.path.append("/home/kalabharath/zinr/main")
import utility.io_util    as io
import utility.ss_util    as ss
import utility.PCSmap     as PCSmap
import utility.ContactMap as contact



def matchSeq2SS(aa_seq, ssfile):
    print aa_seq

    print ssfile
    raw_ss = []
    with open(ssfile) as fin:
        lines = fin.readlines()

    for i in range (0,len(lines)):
        if lines[i] == 'FORMAT %4d %1s %2d %2d %8.3f %8.3f %8.3f %4.2f %s\n':
            print lines[i]
            j = i
            for j in range(j,len(lines)):
                content = lines[j].split()
                if len(content) == 9:
                    raw_ss.append(content)
            break

    print len(aa_seq), len(raw_ss)
    diff = len(raw_ss)-len(aa_seq)

    if diff > 0:
        for i in range(0,diff):
            t_aa=''
            t_ss=''
            for j in range(i,i+len(aa_seq)):
                t_aa=t_aa+raw_ss[j][1]
                t_ss=t_ss+raw_ss[j][-1]
            if t_aa == aa_seq:
                print t_aa, len(t_aa)
                #REMARK     h-Helix    e-Strand   c-Coil (Sequence based)
                t_ss = t_ss.replace('c','L')
                t_ss = t_ss.replace('e','E')
                t_ss = t_ss.replace('h','H')
                print t_ss, len(t_ss)
                return t_ss
    return ss_seq




data = io.readInputDataFiles('input_data.txt')

print data

datatypes = data.keys()
handle, aa_seq = io.readFasta(data['fasta_file'])


ss_seq = matchSeq2SS(aa_seq, data['ss_file'])
#ss_seq = io.readPsiPred(psipred_file)

print ss_seq
ss_def, ss_combi = ss.genSSCombinations(ss_seq)

print ss_combi
io.dumpPickle("ss_profiles.pickle", ss_combi)

# Read in contacts at a given confidence level
if 'contacts_file' in datatypes:
    contacts, contacts_seq = io.readContacts(contactsfile, probability=0.7)


# ss_element format = [ss_type,len_ss,l_loop,r_loop,start,end]
#rank_ss = contact.getContactRoute(ss_def, contacts_seq)

#print rank_ss

# read in PCS data from .npc file from Rosetta's broker file format
pcsdata = io.getPcsTagInfo(ss_seq, data['pcs_broker'])

map_route = PCSmap.getRoute(ss_seq, pcsdata)
print map_route

io.dumpPickle("pcs_route.pickle", map_route)
#io.dumpPickle("contact_route.pickle", rank_ss)

data_dict = {'ss_seq': ss_seq, 'pcs_data': pcsdata, 'aa_seq': aa_seq}


io.dumpPickle("exp_data.pickle", data_dict)
