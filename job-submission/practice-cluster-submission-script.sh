#!/bin/bash 
#$ -S /bin/bash
#$ -V
#$ -N MCNP
#$ -cwd
#$ -j y
#$ -pe orte 8

# Increase the stacksize 
ulimit -s unlimited

export DATAPATH=/share/apps/mcnp/MCNP_DATA
EXE=/share/apps/mcnp/MCNP_CODE/MCNP611/bin/mcnp6.mpi
/opt/openmpi/bin/mpirun $EXE i="$1.inp" o="$1.out" srctp="$1.srctp" mctal="$1.mctal" runtpe="$1.runtpe" > $2