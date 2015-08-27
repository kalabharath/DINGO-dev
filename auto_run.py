import sys, os
sys.path.append("/home/kalabharath/zinr/main")

import utility.io_util as io


run_seq = io.readPickle('pcs_route.pickle')

print run_seq

for i in range(0,len(run_seq)):
    if i == 0:
        mpirun = "mpirun --hostfile mpi_hostfile -np 128 "
        run1 = " python stage1_mpi_run.py | tee stage1.log "
        print mpirun+run1
        os.system(mpirun+run1)
    if i == 1:
        mpirun = "mpirun --hostfile mpi_hostfile -np 128 "
        run1 = " python stage2_mpi_run.py | tee stage2.log "
        print mpirun+run1
        os.system(mpirun+run1)

    if i > 1:
        mpirun = "mpirun --hostfile mpi_hostfile -np 128 "
        run1 = " python stage3x_mpi_run.py | tee stage3_"+str(i)+".log "
        print mpirun+run1
        os.system(mpirun+run1)
