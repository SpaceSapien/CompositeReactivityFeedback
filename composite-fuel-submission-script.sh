#!/bin/bash 
#$ -S /bin/bash
#$ -V
#$ -N MCNP
#$ -cwd
#$ -j y
#$ -pe orte 8

export DATAPATH=/share/apps/mcnp/MCNP_DATA
EXE=/share/apps/mcnp/MCNP_CODE/MCNP6/bin/mcnp6.mpi
mpirun $EXE i=$1 o=$2 > $3