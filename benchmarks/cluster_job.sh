#!/bin/bash
# Test suite and benchmarks on one node (SLURM).  Nodes of partition lambda: 32 cores, 256 GB;
# lambda-hmem (comp1120): 32 cores, 2 TB, for large CCSD/CAS runs, not needed here.
#
# Once, on the login node (Intel's Python; the Fortran compiler: gfortran, or ifx via oneAPI):
#   PY=/opt/intel/oneapi/intelpython/python3.12/bin/python
#   $PY -m venv ~/thistle-env && source ~/thistle-env/bin/activate
#   pip install ".[test]"            # gfortran build (gfortran on the PATH, or a GNU module loaded)
#   # ifx build, in a second environment:
#   #   $PY -m venv ~/thistle-env-ifx && source ~/thistle-env-ifx/bin/activate
#   #   source /opt/intel/oneapi/setvars.sh && FC=ifx pip install ".[test]"
# Then, from the unpacked thistle folder:
#   sbatch benchmarks/cluster_job.sh                # gfortran build -> gnu.npz
#   TAG=ifx THISTLE_ENV=~/thistle-env-ifx sbatch benchmarks/cluster_job.sh   # ifx build -> ifx.npz
#   python benchmarks/cluster_check.py compare gnu.npz ifx.npz          # afterwards: < 1e-8
#SBATCH --job-name=thistle-tests
#SBATCH --partition=lambda          # or lambda-hmem (comp1120)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32          # the whole node: thread scaling up to 32
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=thistle_%j.out
##SBATCH --exclusive               # remove one # for a node to yourself (cleaner timings)

if [ "${TAG:-gnu}" = ifx ]; then
    source /opt/intel/oneapi/setvars.sh > /dev/null
else
    module load gcc 2>/dev/null || true      # if the cluster has no such module, gfortran must be on the PATH
fi
source "${THISTLE_ENV:-$HOME/thistle-env}/bin/activate"
cd "${SLURM_SUBMIT_DIR:-$PWD}"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=1            # Intel Python's NumPy uses MKL; keep it out of the kernel timings

lscpu | head -20
which python && python --version
thistle selftest
PYXSCAT_SLOW=1 python -m pytest -q
python benchmarks/cluster_check.py reference "${TAG:-gnu}.npz"
python benchmarks/cluster_check.py scaling
