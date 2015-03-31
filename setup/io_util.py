#!/usr/bin/env python
import sys, os
"""
Project_Name: setup, File_name: util.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 19/03/15 , Time:3:24 PM
Utility file for frequently used functions
"""
def readFile(filename):
	"""
	read file and return the lines as an array
	:param filename:
	:return:
	"""
	with open(filename,'r') as fin:
		lines = fin.readlines()
	return  lines

def writeFile(array_of_lines, filename):
	"""
	write output to a file
	:param array_of_lines:
	:param filename:
	:return:
	"""
	fout = open(filename,'w')
	for line in array_of_lines:
		fout.write(line)
		fout.write("\n")
	fout.close()
def readFasta(filename):
	"""
	reads in FastA file and returns seq output
	:param filename:
	:return:
	"""
	with open (filename,'r') as fin:
		lines=fin.readlines()

	header = lines[0].rstrip()
	seq = ''
	for i in range(1,len(lines)):
		seq=seq+lines[i].rstrip()
	return header, seq

def readPsiPred(filename):
	"""
	read in ss_seq as predicted by PsiPred
	:param filename:
	:return: ss_seq
	"""
	with open(filename,'r') as fin:
		lines = fin.readlines()

	ss_seq=''

	for line in lines:
		if line[0:5]=='Pred:':
			ss_seq = ss_seq+line[6:-1]
	return  ss_seq

def readContacts(filename, probability):
	"""
    Parses MetaPSICOV's contact file and returns the
	:param filename:
	:param probability:
	:return contacts array
	"""

	with open(filename,'r') as fin:
		lines = fin.readlines()
	contacts = []
	contacts_seq={}
	for line in lines:
		if line != '\n':
			res1, res2, zero, distance, precision = line.split()
			precision = round(float(precision),1)
			if precision >= probability:
				contacts.append([int(res1), int(res2), float(distance), float(precision)])
				contacts_seq.setdefault(int(res1),[]).append(int(res2))
				contacts_seq.setdefault(int(res2),[]).append(int(res1))
	return contacts, contacts_seq