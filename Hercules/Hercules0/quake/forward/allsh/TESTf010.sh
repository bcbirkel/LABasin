
#!/bin/bash
#BSUB -oo TESTf001.out
#BSUB -eo TESTf001.err
#BSUB -q q_64p_6h  -n 64
# New flags for new planes (04/08)
export MPICH_PTL_SEND_CREDITS=-1
export MPICH_PTL_UNEX_EVENTS=20000
CVM_DBNAME=XXXX
MY_INPUT=/tmpu/lramg_g/lramg/MEXICO/GREEN_DATABASE/tests/f010/input
  
/opt/SC/intel/impi/4.1.0.024/intel64/bin/mpirun -np 64 /tmpu/lramg_g/lramg/MEXICO/GREEN_DATABASE/tests/hercules/quake/forward/psolve \
   $CVM_DBNAME \
   $MY_INPUT/input.in \
   $MY_INPUT/input.in \
   $MY_INPUT/mesh.e \
   $MY_INPUT/disp.out
