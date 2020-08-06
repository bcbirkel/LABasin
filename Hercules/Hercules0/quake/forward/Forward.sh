
#!/bin/bash
#  New flags for new planes (04/08)
% Run setup file
. /usr/usc/openmpi/1.8.7/setup.sh


CVM_DBNAME=XXXX
MY_INPUT=/home/scec-00/birkel/S20140329/Simulation/Forward/Source1/input

mpirun -n 10 /home/scec-00/birkel/Hercules0/quake/forward/psolve \
   $CVM_DBNAME \
   $MY_INPUT/input.in \
   $MY_INPUT/input.in \
   $MY_INPUT/mesh.e \
   $MY_INPUT/disp.out
