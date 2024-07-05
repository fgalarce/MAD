#!/bin/bash
#SBATCH -N 4
 
cpus_per_node=${SLURM_JOB_CPUS_PER_NODE%\(*}

module load gcc/9.3.0

export MAD_ROOT=/home/bebbiofl/mad/source/
export MDATA=/home/bebbiofl/mad/data/mesh/
export PETSC_ROOT=/home/bebbiofl/mad/petsc/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PETSC_ROOT/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PETSC_ROOT/build/lib/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PETSC_ROOT/arch-linux2-c-opt/externalpackages/mpich-3.3.1/lib/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MAD_ROOT/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MAD_ROOT/fem/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MAD_ROOT/simplex/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MAD_ROOT/tools/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MAD_ROOT/io/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MAD_ROOT/cfd/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MAD_ROOT/la/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MAD_ROOT/assimilation/
export LD_LIBRARY_PATH

##cd ../../../
##make
##cd projects/pdes/PEDup_rectangle/
make
/home/bebbiofl/mad/petsc/build/bin/mpirun -np 100 ./poro.exe par_toy
