#!/bin/bash
#SBATCH -o job.%j.out
#SBATCH -p 640
#SBATCH -J C
#SBATCH --nodes=1
#SBATCH -x cpu10
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
python ../basis_opt/main.py
#abacus_work=/home/ghj/abacus/abacus_abfs/abacus-develop/build/abacus
#IB='-env OMP_NUM_THREADS=48'
#mpirun -n 1 -ppn 1 $IB $abacus_work >  Ne.out
#mpirun -n 1 -ppn 1 $IB /home/ghj/abacus/LibRPA/build/chi0_main.exe 16 0 > LibRPA_Ne.out

echo End Time: `date`
