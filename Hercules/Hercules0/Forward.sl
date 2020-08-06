#!/bin/bash -l

# Email messages
#SBATCH --mail-user=birkel@usc.edu
#SBATCH --mail-type=FAIL

% Run setup file
. /usr/usc/openmpi/1.8.7/setup.sh

# Run command
MY_INPUT=/home/scec-00/birkel/${E}/Simulation/${F}/${R}/input
mpirun -n $T /home/scec-00/birkel/Hercules${H}/quake/forward/psolve XXXX \
$MY_INPUT/input.in \
$MY_INPUT/input.in \
$MY_INPUT/mesh.e \
$MY_INPUT/disp.out
