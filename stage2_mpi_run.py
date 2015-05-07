#!/usr/bin/env python

"""
Project_Name: main, File_name: stage_1_mpi_run.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 13/04/15 , Time:10:05 AM

Perform stage 1 in perfect parallel
"""
#TODO write detailed comments

from   mpi4py import  MPI
import utility.stage2_util as util


#Define MPI messaage tags
tags = util.enum('READY', 'DONE', 'EXIT', 'START')

#Init MPI comm

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()
status = MPI.Status()

##

if rank == 0:
    tasks = util.getRunSeq()