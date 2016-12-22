#!/bin/bash 
#$ -S /bin/bash
#$ -V
#$ -N MCNP
#$ -cwd
#$ -j y
#$ -pe orte 72

# DISPLAY for MCNP executables
DISPLAY=":0.0"

# Increase the stacksize 
ulimit -s unlimited

LD_LIBRARY_PATH="/share/apps/gcc/gcc-6.2.0/lib64"
DATAPATH=/share/apps/mcnp/MCNP/MCNP_DATA
EXE=/share/apps/mcnp/MCNP/MCNP_CODE/MCNP611/bin/mcnp6.mpi
MPI=/opt/openmpi/bin/mpirun

$MPI $EXE i="$1.inp" o="$1.out" srctp="$1.srctp" mctal="$1.mctal" runtpe="$1.runtpe" > $2
