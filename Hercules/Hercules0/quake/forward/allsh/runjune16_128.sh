#!/bin/bash
#BSUB -oo june16.out
#BSUB -eo june15.err
#BSUB -q q_128p_12h -n 128
# New flags for new planes (04/08)
export MPICH_PTL_SEND_CREDITS=-1
export MPICH_PTL_UNEX_EVENTS=20000
CVM_DBNAME=XXXX

MY_INPUT=/tmpu/lramg_g/lramg/MEXICO/RUNS/June16thEartqhake/input/
 
/opt/SC/intel/impi/4.1.0.024/intel64/bin/mpirun -np 128 /tmpu/lramg_g/lramg/MEXICO/GREEN_DATABASE/tests/hercules/quake/forward/psolve \
   $CVM_DBNAME \
   $MY_INPUT/input.in \
   $MY_INPUT/input.in \
   $MY_INPUT/mesh.e \
   $MY_INPUT/disp.out
