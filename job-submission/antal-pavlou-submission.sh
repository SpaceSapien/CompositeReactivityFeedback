#!/bin/bash 
#$ -S /bin/bash
#$ -V
#$ -N MCNP
#$ -cwd
#$ -j y
#$ -pe orte 64

# Increase the stacksize 
ulimit -s unlimited

export DATAPATH=/share/apps/MCNP-ANDREW-PAVLOU/MCNP_DATA
export DATAPATH
EXE=/share/apps/MCNP-ANDREW-PAVLOU/MCNP_CODE/MCNP6/bin/mcnp6.mpi
MPI="/opt/openmpi/bin/mpirun"
$MPI $EXE i="$1.inp" o="$1.out" srctp="$1.srctp" runtpe="$1.runtpe" mctal="$1.mctal" > $2