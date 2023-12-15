#!/bin/bash
#SBATCH -o job.%j.out
#SBATCH -p 640
#SBATCH -J test
#SBATCH --nodes=1
##SBATCH -x cpu01
##SBATCH --ntasks=16
#SBATCH --cpus-per-task=48
#SBATCH --ntasks-per-node=1

#module load gcc10.2
#module load intel20u4

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo Working directory is $PWD
echo This job has allocated $SLURM_JOB_CPUS_PER_NODE cpu cores.

echo Begin Time: `date`
workdir=$(basename `pwd`)
### * * * Running the tasks * * * ###
abacus_work=/home/ghj/abacus/230820/abacus-develop/build/abacus
IB='-env OMP_NUM_THREADS=48'
mpirun -n 1 -ppn 1 $IB $abacus_work > $workdir.$SLURM_JOB_ID.out
#mpirun /home/ghj/abacus/LibRPA/build/chi0_main.exe 14 0 > LibRPA_$workdir.$SLURM_JOB_ID.out

echo End Time: `date`
