#!/usr/bin/env python

"""

Project_Name: main, File_name: master_mpi
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 13/04/15 , Time:10:05 AM

One master mpi and one master_search to eliminate redundancy across four files
should be called as a part of the sequence of smotif assembly files

"""

# sys.path.append('../../main/')
import argparse
import time
import traceback
from   mpi4py import MPI

import ranking.SmotifRanking as srank
import smotif_search as msearch
import utility.masterutil as mutil
import utility.stage2_util as util

# Define MPI message tags

tags = mutil.enum('READY', 'DONE', 'EXIT', 'START')
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()
status = MPI.Status()


def killall(processes):
    """
    Kill all the subprocess when requested
    :param processes:
    :return: True or False
    """
    for i in range(0, processes - 1):
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        source = status.Get_source()
        tag = status.Get_tag()
        if tag == tags.READY:
            comm.send(None, dest=source, tag=tags.EXIT)
    return True


#####################################  Define cmd line argument parser #############################################

parser = argparse.ArgumentParser(description='DINGO-PCS Master MPI process that manages all jobs.')
parser.add_argument('--stage', type=int, help='specify the stage of  the Smotif assembly')
parser.add_argument('--numhits', type=int, help='Top number of hits to be selected from previous assembly')
args = parser.parse_args()
#####################################  Define cmd line argument parser #############################################

# Rank '0' specifies the master process

if rank == 0:

    ##################################  Extract top hits from previous stage ########################################

    sse_index = 9999999

    if args.stage == 1:
        tasks = mutil.getRunSeq()  # there are no hits to extract if it is the 1st stage

    else:
        # for 2nd stage and beyond
        try:
            try:
                # Restart from the already assembled top hits if possible
                tasks, sse_index = util.start_top_hits(args.numhits, args.stage)
            except:
                # Assemble top hits from the previously generated hits
                tasks, sse_index = srank.getRunSeq(args.numhits, args.stage)
        except:
            # print what went wrong and terminate the slave processes
            traceback.print_exc()
            print "Couldn't extract top hits within the specified cutoffs: Exiting..."
            killall(size)
            exit()

        if sse_index == 9999999:
            # kill all slaves if there is there is EOL
            # only makes sense for self submitting jobs
            killall(size)
            exit()

    ##################################  Generate and distribute job index array ########################################

    stime = time.time()

    # print tasks, len(tasks) # this will be the new tasks
    task_index = 0  # control the number of processes with this index number
    finished_task = 0
    num_workers = size - 1  # 1 processor is reserved for master.
    closed_workers = 0  # control the workers with no more work that can be assigned

    print ("Master starting with {} workers".format(num_workers))
    while closed_workers < num_workers:
        # Manage/distribute all processes in this while loop
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        source = status.Get_source()
        tag = status.Get_tag()
        if tag == tags.READY:
            # worker process is ready, send some task to do.
            if task_index < len(tasks):
                comm.send([tasks[task_index], args.stage], dest=source, tag=tags.START)
                task_index += 1  # increment its
            else:
                # everything is done, send exit signal
                comm.send(None, dest=source, tag=tags.EXIT)
        elif tag == tags.DONE:
            # take the result from the worker
            if data:
                print "Receiving data"
            ctime = time.time()
            elapsed = ctime - stime
            finished_task += 1
            print "Finishing..", finished_task, "of", len(tasks), "Smotifs, Elapsed", round((elapsed) / (60), 2), "mins"
        elif tag == tags.EXIT:
            closed_workers += 1
    # Rename temprary files
    util.rename_pickle(sse_index)
    print "All Done, Master exiting"
    exit()


# On the worker processes
else:

    while True:  # initiaite infinite loop
        comm.send(None, dest=0, tag=tags.READY)
        # Signal the master process that you are READY

        task = comm.recv(source=0, tag=MPI.ANY_SOURCE, status=status)
        tag = status.Get_tag()
        if tag == tags.START:
            result = False
            if args.stage == 1:
                result = msearch.S1SmotifSearch(task)
            else:
                result = msearch.sXSmotifSearch(task)

            comm.send(result, dest=0, tag=tags.DONE)
        elif tag == tags.EXIT:
            # break the infinite loop because there is no more work that can be assigned
            break

    # Signal EXIT to the master process
    comm.send(None, dest=0, tag=tags.EXIT)
