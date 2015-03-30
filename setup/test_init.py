#!/usr/bin/env python
import sys, os
import  io_util as io_ut
import ss_util as ss_ut
"""
Project_Name: setup, File_name: test_init.py 
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 19/03/15 , Time:4:23 PM
"""

header, seq = io_ut.readFasta("seq.fasta")
ss_seq = io_ut.readPsiPred('2m47.psipass2')
contacts = io_ut.readContacts('./2m47_metapsicov_contacts.txt', 0.9)
print seq, len(seq)
print ss_seq, len(ss_seq)
# print contacts, len(contacts)
ss_def, ss_combi = ss_ut.genSSCombinations(ss_seq)

print ss_def
print ss_combi
