#!/usr/bin/env python

"""
Project_Name: main, File_name: stage_0.py 
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 31/03/15 , Time:10:56 AM

Prepare all the relavant files for stage1 & 2
"""

from setup.io_util import *
from setup.ss_util import *
from setup.PCSmap import *
from setup.ContactMap import *


fasta_file = './2m47.fasta'
psipred_file = './2m47.psipass2'
contactsfile = './2m47_metapsicov_contacts.txt'
pcs_broker = './broker-ts3.txt'



ss_seq = readPsiPred(psipred_file)
print ss_seq
ss_def, ss_combi = genSSCombinations(ss_seq)
print ss_def
contacts, contacts_seq = readContacts(contactsfile,0.7)
print contacts_seq
#ss_element format = [ss_type,len_ss,l_loop,r_loop,start,end]


rank_ss = getContactRoute(ss_def, contacts_seq)

"""
pcsdata = getPcsTagInfo(ss_seq,pcs_broker)
print pcsdata
map_route = getRoute(ss_seq,pcsdata)

print map_rout
"""