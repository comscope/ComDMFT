#!/bin/sh
#SBATCH -p long
#SBATCH -t 06:00:00
#SBATCH -A ms17q1
#SBATCH -N 4
#SBATCH --ntasks-per-node=36
#SBATCH -o hostname_%j.out
#SBATCH -e hostname_%j.err

cd $SLURM_SUBMIT_DIR

module load intel/PSXE2017.u4

export OMP_NUM_THREADS=1

echo begin at: `date` >> logfile

srun -n 144 ~/V101/ctqmc/CPU/CTQMC params >& screen.out

echo end at: `date` >> logfile
