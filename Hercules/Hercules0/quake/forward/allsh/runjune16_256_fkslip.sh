#!/bin/bash
#BSUB -oo june16_fkslip.out
#BSUB -eo june16_fkslip.err
#BSUB -q q_256p_12h -n 256
# New flags for new planes (04/08)
export MPICH_PTL_SEND_CREDITS=-1
export MPICH_PTL_UNEX_EVENTS=20000
CVM_DBNAME=XXXX

MY_INPUT=/tmpu/lramg_g/lramg/MEXICO/RUNS/June16thEartqhake_fakeslip/input/
 
/opt/SC/intel/impi/4.1.0.024/intel64/bin/mpirun -np 256 /tmpu/lramg_g/lramg/MEXICO/GREEN_DATABASE/tests/hercules/quake/forward/psolve \
   $CVM_DBNAME \
   $MY_INPUT/input.in \
   $MY_INPUT/input.in \
   $MY_INPUT/mesh.e \
   $MY_INPUT/disp.out
