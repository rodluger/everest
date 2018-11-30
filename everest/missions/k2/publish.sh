#!/bin/sh

# Prevent numpy from multithreading
export OPENBLAS_NUM_THREADS=1

# CD into our working directory
cd ${EVERESTDAT}

if [ $NODES -eq 1 ]; then
  # Run on a single node with Python multiprocessing
  python -c "import everest; everest.missions.k2.slurm._Publish($CAMPAIGN, $SUBCAMPAIGN, $STRKWARGS)"
else
  # Run on multiple nodes with MPI
  mpirun -np $SLURM_NTASKS python -c "import everest; everest.missions.k2.slurm._Publish($CAMPAIGN, $SUBCAMPAIGN, $STRKWARGS)"
fi