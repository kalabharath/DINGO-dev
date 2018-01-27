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
from mpi4py import MPI
from ranking.NoeStageRank import *
import alt_smotif_search as alt_search
import utility.masterutil as mutil
import utility.io_util as io
import utility.alt_smotif_util as altutil

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
    count = 0
    while True:
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        source = status.Get_source()
        tag = status.Get_tag()
        if tag == tags.READY:
            comm.send(None, dest=source, tag=tags.EXIT)
            count += 1
        if count == processes -1:
            break
    return True


# ********************* Define cmd line argument parser *********************

parser = argparse.ArgumentParser(description='DINGO-Refine Master MPI process that manages all jobs.')
parser.add_argument('--infile', type=int, help='specify the top_hits file')
parser.add_argument('--stage', type=int, help='specify the stage of  the Smotif assembly')
parser.add_argument('--numhits', type=int, help='Top number of hits to be selected')
args = parser.parse_args()

# *********************   Define cmd line argument parser *********************

# Rank '0' specifies the master process

if rank == 0:

    # *********************   Extract top hits *********************

    in_file = str(args.infile)+"_tophits.gzip"
    print "infile ", in_file

    try:
        tasks = io.readGzipPickle(in_file)
        print "len of tasks", len(tasks)
    except:
        traceback.print_exc()
        print "There are no entries in the tophits file, nothing to refine"
        killall(size)
        exit()

    # ********************* Generate and distribute job index array *********************

    stime = time.time()

    try:
        if len(tasks):
            pass
    except:
        print "killing all processes!"
        killall(size)
        exit()

    # print tasks, len(tasks) # this will be the new tasks
    task_index = 0  # control the number of processes with this index number
    finished_task = 0
    num_workers = size - 1  # 1 processor is reserved for master.
    closed_workers = 0  # control the workers with no more work that can be assigned

    print ("Master starting with {} workers".format(num_workers))
    total_data = []
    for entry in tasks:
        total_data.append(entry)
    lowest_noe_energy = altutil.get_lowest_noe_energy(tasks)
    print "Average lowest NOE energy is :", lowest_noe_energy

    total_jobs, alt_sse_profile = altutil.compute_jobs(tasks)

    while closed_workers < num_workers:

        # Manage/distribute all processes in this while loop
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        source = status.Get_source()
        tag = status.Get_tag()
        if tag == tags.READY:
            # worker process is ready, send some task to do.
            if task_index < len(total_jobs):
                t_job = total_jobs[task_index]
                send_job = [tasks[t_job[0]], alt_sse_profile[t_job[1]], args.stage, task_index, lowest_noe_energy]
                comm.send(send_job, dest=source, tag=tags.START)
                task_index += 1  # increment its
            else:
                # everything is done, send exit signal
                comm.send(None, dest=source, tag=tags.EXIT)
        elif tag == tags.DONE:
            # take the result from the worker
            if data:
                for hit in data:
                    total_data.append(hit)
            ctime = time.time()
            elapsed = ctime - stime
            finished_task += 1
            print "Finishing..", finished_task, "of", len(total_jobs), "Smotifs, Elapsed", round((elapsed) / (60), 2), "mins"
        elif tag == tags.EXIT:
            closed_workers += 1

    # consolidate top_hits and dump files here
    print "Total number of hits  found are : ",len(total_data)
    # ranked_data = rank_assembly(total_data, args.numhits)

    ranked_data = rank_assembly_with_clustering(total_data, args.numhits)
    print len(ranked_data)
    io.dumpGzipPickle(str(args.infile) + "_refined_tophits.gzip", ranked_data)
    print "All Done, Master exiting"
    exit()


# On the worker processes
else:

    while True:  # initiate infinite loop
        comm.send(None, dest=0, tag=tags.READY)
        # Signal the master process that you are READY
        work = comm.recv(source=0, tag=MPI.ANY_SOURCE, status=status)
        tag = status.Get_tag()
        if tag == tags.START:
            result = alt_search.altSmotifSearch(work)
            comm.send(result, dest=0, tag=tags.DONE)
        elif tag == tags.EXIT:
            # break the infinite loop because there is no more work that can be assigned
            break

    # Signal EXIT to the master process
    comm.send(None, dest=0, tag=tags.EXIT)
