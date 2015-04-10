#!/usr/bin/env python

"""
Project_Name: main, File_name: stage_0.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 31/03/15 , Time:10:56 AM

Prepare all the relavant files for stage1 & 2


from setup.io_util import *
from setup.ss_util import *
from setup.PCSmap import *
from setup.ContactMap import *
"""

import setup.io_util as io
import  setup.ss_util as ss
import  setup.PCSmap as PCSmap
import  setup.ContactMap as contact

fasta_file = './2m47.fasta'
psipred_file = './2m47.psipass2'
contactsfile = './2m47_metapsicov_contacts.txt'
pcs_broker = './broker-ts3.txt'



ss_seq = io.readPsiPred(psipred_file)
print ss_seq
ss_def, ss_combi = ss.genSSCombinations(ss_seq)
print ss_def
contacts, contacts_seq = io.readContacts(contactsfile,0.7)
print contacts_seq
#ss_element format = [ss_type,len_ss,l_loop,r_loop,start,end]


rank_ss = contact.getContactRoute(ss_def, contacts_seq)

print rank_ss
"""
pcsdata = getPcsTagInfo(ss_seq,pcs_broker)
print pcsdata
map_route = getRoute(ss_seq,pcsdata)

print map_rout
"""
