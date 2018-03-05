#!/usr/bin/env python

"""

Project_Name: DINGO, File_name: gather_and_stitch
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 4/03/18 , Time:11:00 AM

"""
import  os, glob, argparse
import utility.io_util as io


def get_top_hits_file():
    start = 1
    while True:
        tfile = str(start)+"_refined_tophits.gzip"
        if os.path.isfile(tfile):
            start +=1
        else:
            return tfile, start


def combine_data():
    hits = []

    file_list = glob.glob("rtx_*.gzip")
    if file_list:
        pass
    else:
        return False
    for f in file_list:
        print "Extracting data from file..", f
        t_hits = io.readGzipPickle(f)
        for t_hit in t_hits:
            hits.append(t_hit)
    return hits


def gather_and_stitch(seq, tfile):

    total_data =[]
    in_file  = str(seq)+"_tophits.gzip"

    if os.path.isfile(in_file):
        print "Extracting data from non_refined file:", in_file
        tasks = io.readGzipPickle(in_file)
        for entry in tasks:
            total_data.append(entry)
    else:
        print "Something is wrong in the datafile: ",in_file

    refined_data = combine_data()
    if refined_data:
        for entry in refined_data:
            total_data.append(entry)
    else:
        pass

    from ranking.NoeStageRank import rank_assembly_with_clustering
    ranked_data = rank_assembly_with_clustering(total_data, args.numhits)
    io.dumpGzipPickle(str(tfile), ranked_data)
    # delete files

    try:
        print "Deleting old rtx_* files!"
        rm_files = "rm rtx_*.gzip"
        os.system(rm_files)
    except NameError:
        print "No rtx_* files exist!"
    return True


if __name__ == "__main__":

    # ********************* Define cmd line argument parser *********************

    parser = argparse.ArgumentParser(description='Gather and stitch to recover from Memory overflows.')
    parser.add_argument('--infile', type=int, help='specify the top_hits file')
    parser.add_argument('--stage', type=int, help='specify the stage of  the Smotif assembly')
    parser.add_argument('--numhits', type=int, help='Top number of hits to be selected')
    args = parser.parse_args()

    # *********************   Define cmd line argument parser *********************
    print "Gather and stitching fragmented data! This May take a long time \n"
    seq = args.infile
    tfile = str(args.infile)+"_refined_tophits.gzip"
    #tfile, seq = get_top_hits_file()
    gather_and_stitch(seq, tfile)
