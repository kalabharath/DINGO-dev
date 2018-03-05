import sys, os, glob

run1 = "#!/bin/bash \n #PBS -P xc4 \n #PBS -q normal \n #PBS -l walltime=5:00:00 \n #PBS -l mem=256GB \n #PBS -l ncpus=128 \n #PBS -l wd \n"

run2 = "#!/bin/bash \n #PBS -P xc4 \n #PBS -q normal \n #PBS -l walltime=5:00:00 \n #PBS -l mem=256GB \n #PBS -l ncpus=128 \n #PBS -l wd \n"

run3 = "#!/bin/bash \n #PBS -P xc4 \n #PBS -q normal \n #PBS -l walltime=5:00:00 \n #PBS -l mem=32GB \n #PBS -l ncpus=1 \n #PBS -l wd \n"

with open("run.sh") as fin:
    lines = fin.readlines()

count = 0
for line in lines():
    if "alt" in line:
        count += 1


for i in range(count):

    if i == 0:
        file_name = "sub_" + str(i) + ".sh"
        fout = open(file_name, 'w')
        fout.write(run1)
        job1 = "mpirun -np 128 python ../../main/Dingo_mpi.py --stage 1  --numhits 127 \n"
        fout.write(job1)
        job2 = "python inter_rmsd " + str(i) + " >" + str(i) + ".log \n"
        fout.write(job2)

    if i == 1:
        job1 = "mpirun -np 128 python ../../main/Dingo_mpi.py --stage 2  --numhits 127 \n"
        fout.write(job1)
        job2 = "python inter_rmsd " + str(i) + " >" + str(i) + ".log \n"
        fout.write(job2)
        job3 = "mpirun - np 128 python.. /..a/ main / Dingo_alt_mpi.py - -infile 1 - -stage 2 - -numhits 127\n"
        fout.write(job3)
        fgather = "gsub_" + str(i) + ".sh\n"
        qsub = "qsub " + fgather
        fout.write(qsub)
        fout.close()


        fout = open(fgather, 'w')
        fout.write(run3)
        jobf = "python ../../main/gather_and_stitch.py --infile" + str(i) + " --numhits 127\n"
        fout.write(jobf)
        job2 = "python inter_rmsd " + str(i) + " >" + str(i) + ".refined_log \n"
        fout.write(job2)
        tfile_name = "sub_" + str(i + 1) + ".sh\n"
        qsub = "qsub " + tfile_name
        fout.write(qsub)
        fout.close()

    if i > 1:
        file_name = "sub_" + str(i) + ".sh"
        fout = open(file_name, 'w')
        fout.write(run1)
        job1 = "mpirun -np 128 python ../../main/Dingo_mpi.py --stage 3  --numhits 127 \n"
        fout.write(job1)
        job2 = "python inter_rmsd " + str(i) + " >" + str(i) + ".log \n"
        fout.write(job2)
        job3 = "mpirun - np 128 python.. /..a/ main / Dingo_alt_mpi.py - -infile"+str(i)+" - -stage 2 - -numhits 127\n"
        fout.write(job3)
        fgather = "gsub_" + str(i) + ".sh\n"
        qsub = "qsub " + fgather
        fout.close()

        fout = open(fgather, 'w')
        fout.write(run3)
        jobf = "python ../../main/gather_and_stitch.py --infile" + str(i) + " --numhits 127\n"
        fout.write(jobf)
        job2 = "python inter_rmsd " + str(i) + " >" + str(i) + ".refined_log \n"
        fout.write(job2)
        tfile_name = "sub_" + str(i + 1) + ".sh\n"
        qsub = "qsub " + tfile_name
        fout.write(qsub)
        fout.close()