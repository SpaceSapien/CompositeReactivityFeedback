#!/bin/bash 
#$ -S /bin/bash
#$ -V
#$ -N MCNP
#$ -cwd
#$ -j y
#$ -pe orte 36

# Increase the stacksize 
ulimit -s unlimited

DATAPATH="/share/apps/mcnp/MCNP_DATA"
export DATAPATH
EXE="/share/apps/mcnp/MCNP_CODE/MCNP611/bin/mcnp6.mpi"
MPI="/opt/openmpi/bin/mpirun"
LD_LIBRARY_PATH="/share/apps/gcc-5-3-0/lib64:$LD_LIBRARY_PATH"
$MPI $EXE i="$1.inp" o="$1.out" srctp="$1.srctp" mctal="$1.mctal" runtpe="$1.runtpe" > $2