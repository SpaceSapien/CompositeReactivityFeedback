#!/bin/bash 
#$ -S /bin/bash
#$ -V
#$ -N MCNP
#$ -cwd
#$ -j y
#$ -pe orte 8

DATAPATH="/share/apps/mcnp/MCNP_DATA"
EXE="/share/apps/mcnp/MCNP_CODE/MCNP6/bin/mcnp6.mpi"
MPI="/share/apps/openmpi-1.10.2/bin/mpirun"
LD_LIBRARY_PATH="/share/apps/gcc-5-3-0/lib64:$LD_LIBRARY_PATH"
MCNP_INPUT="i=\"$1.inp\" o=\"$1.out\" srctp=\"$1.srctp\" mctal=\"$1.mctal\" runtpe=\"$1.runtpe\" > $2"
COMMAND="$MPI $EXE $MCNP"

echo $COMMAND

#run the command
$($COMMAND)
