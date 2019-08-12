#!/bin/sh

# CD into our working directory
cd ${EVERESTDAT}

# Run on a single node with Python multiprocessing
python -c "import everest; everest.missions.k2.slurm._Download($CAMPAIGN, $SUBCAMPAIGN)"