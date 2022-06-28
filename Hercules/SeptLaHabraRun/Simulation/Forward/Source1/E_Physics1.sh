#!/bin/bash

module load intel/19.0.4 intel-mpi

MY_INPUT=/project/scec_608/luisalbe/S20200604/Simulation/Forward/Source1/input
srun --mpi=pmi2 -n 8  /home1/luisalbe/Hercules0/quake/forward/psolve XXXX \
$MY_INPUT/input.in \
$MY_INPUT/input.in \
$MY_INPUT/mesh.e \
$MY_INPUT/disp.out
