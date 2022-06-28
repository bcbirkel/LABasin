

#!/bin/bash
#BSUB -oo st.out
#BSUB -eo st.err
#BSUB -q q_512p_24h -n 512
# New flags for new planes (04/08)
export MPICH_PTL_SEND_CREDITS=-1
export MPICH_PTL_UNEX_EVENTS=20000
CVM_DBNAME=XXXX

for i in 51 52 53 54
do

MY_INPUT=/tmpu/lramg_g/lramg/MEXICO/GREEN_DATABASE/IISSNstations/station_$i/001/input/

 
/opt/SC/intel/impi/4.1.0.024/intel64/bin/mpirun -np 512 /tmpu/lramg_g/lramg/MEXICO/GREEN_DATABASE/tests/hercules/quake/forward/psolve \
   $CVM_DBNAME \
   $MY_INPUT/input.in \
   $MY_INPUT/input.in \
   $MY_INPUT/mesh.e \
   $MY_INPUT/disp.out

MY_INPUT=/tmpu/lramg_g/lramg/MEXICO/GREEN_DATABASE/IISSNstations/station_$i/010/input/
  
/opt/SC/intel/impi/4.1.0.024/intel64/bin/mpirun -np 512 /tmpu/lramg_g/lramg/MEXICO/GREEN_DATABASE/tests/hercules/quake/forward/psolve \
   $CVM_DBNAME \
   $MY_INPUT/input.in \
   $MY_INPUT/input.in \
   $MY_INPUT/mesh.e \
   $MY_INPUT/disp.out

MY_INPUT=/tmpu/lramg_g/lramg/MEXICO/GREEN_DATABASE/IISSNstations/station_$i/100/input/
  
/opt/SC/intel/impi/4.1.0.024/intel64/bin/mpirun -np 512 /tmpu/lramg_g/lramg/MEXICO/GREEN_DATABASE/tests/hercules/quake/forward/psolve \
   $CVM_DBNAME \
   $MY_INPUT/input.in \
   $MY_INPUT/input.in \
   $MY_INPUT/mesh.e \
   $MY_INPUT/disp.out

done