#!/bin/bash
#SBATCH --job-name=mrci-water
#SBATCH --partition=lambda
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=mrci_water_%j.out

export MOLCAS=/home/lambda-data/archive/newc6739/Programs/OpenMolcas/OpenMolcas_v24_10_2024_11_15_gnu_mkl/
export PATH=$MOLCAS:$MOLCAS/bin:$PATH
export MOLCAS_MEM=8000                   # MB
export MOLCAS_NPROCS=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MOLCAS_WORKDIR=$SLURM_SUBMIT_DIR/work   # scratch next to the input, so the files stay
export MOLCAS_KEEP_WORKDIR=YES
export MOLCAS_PROJECT=water_mrci
mkdir -p "$MOLCAS_WORKDIR"
cd "$SLURM_SUBMIT_DIR"

pymolcas water_mrci.input > water_mrci.log 2>&1
echo "exit status $?"
ls -la . work/water_mrci 2>/dev/null
