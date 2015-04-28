#!/usr/bin/env python

"""
Project_Name: main, File_name: stage_0.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 31/03/15 , Time:10:56 AM

Prepare all the relavant files for stage1 & 2
"""

import  utility.io_util    as io
import  utility.ss_util    as ss
import  utility.PCSmap     as PCSmap
import  utility.ContactMap as contact

data = io.readInputDataFiles('input_data.txt', ['fasta_file'])

print data


#TODO refactor of reading input files to readInputDataFiles


fasta_file = './2m47.fasta'
psipred_file = './2m47.psipass2'
contactsfile = './2m47_metapsicov_contacts.txt'
pcs_broker = './broker-ts3.txt'

ss_seq = io.readPsiPred(psipred_file)
print ss_seq
ss_def, ss_combi = ss.genSSCombinations(ss_seq)

io.dumpPickle("ss_profiles.pickle", ss_combi)

# Read in contacts at a given confidence level
contacts, contacts_seq = io.readContacts(contactsfile, probability=0.7)


#ss_element format = [ss_type,len_ss,l_loop,r_loop,start,end]
rank_ss = contact.getContactRoute(ss_def, contacts_seq)

print rank_ss

# read in PCS data from .npc file from Rosetta's broker file format
pcsdata = io.getPcsTagInfo(ss_seq,pcs_broker)

map_route = PCSmap.getRoute(ss_seq,pcsdata)
print map_route


io.dumpPickle("pcs_route.pickle", map_route)
io.dumpPickle("contact_route.pickle", rank_ss)
