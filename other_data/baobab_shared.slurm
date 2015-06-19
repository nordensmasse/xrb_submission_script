#!/bin/bash
#SBATCH --mail-user=mads.sorensen@unige.ch
#SBATCH --job-name=DummyJobName
#SBATCH --output=DummyJobName.out
#SBATCH --error=DummyJobName.err
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=8000
#SBATCH --time=12:00:00

export PATH="/home/sorensen/anaconda/bin:$PATH"
export OMP_DYNAMIC=true
export OMP_NUM_THREADS=16

export MESASDK_ROOT=~/mesasdk
source $MESASDK_ROOT/bin/mesasdk_init.sh
export MESA_DIR=~/mesa-r7184

python StartJobRLO.py
