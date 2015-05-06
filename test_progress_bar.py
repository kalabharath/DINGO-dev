#!/usr/bin/env python

"""
Project_Name: main, File_name: stage_1_mpi_run.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 13/04/15 , Time:10:05 AM

Perform stage 1 in perfect parallel
"""
#TODO write detailed comments

from   mpi4py import  MPI
import utility.stage1_util as util
from progress.bar import ChargingBar
import stage_1_search as S1search

#Define MPI messaage tags

tags = util.enum('READY', 'DONE', 'EXIT', 'START')

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()
status = MPI.Status()





# on the master process
if rank == 0:

    tasks = util.getRunSeq()

    #tasks = [[10,0]]
    #tasks = tasks[80:85]

    #print tasks, len(tasks) # this will be the new tasks
    task_index =0 # control the number of processes with this index number
    num_workers = size -1 # 1 processor is reserved for master.
    closed_workers = 0 # control the workers with no more work that can be assigned

    #print ("Master starting with {} workers".format(num_workers))

    ###Progress bar
    import time
    bar = ChargingBar('Processing', max=20)

    for i in range(20):
        bar.next()
        time.sleep(1)
    bar.finish()
