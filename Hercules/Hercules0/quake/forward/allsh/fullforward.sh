
#!/bin/bash
#BSUB -oo proc_fw.out
#BSUB -eo proc_fw.err
#BSUB -q q_512p_6h  -n 512
# New flags for new planes (04/08)
export MPICH_PTL_SEND_CREDITS=-1
export MPICH_PTL_UNEX_EVENTS=20000
CVM_DBNAME=XXXX
MY_INPUT=/tmpu/lramg_g/lramg/CHINA/geodeticfullforward/input/
  
/opt/SC/intel/impi/4.1.0.024/intel64/bin/mpirun -np 512 /tmpu/lramg_g/lramg/CHINA/HERCULES_OLLIN_MAY2014/quake/forward/psolve \
   $CVM_DBNAME \
   $MY_INPUT/input.in $MY_INPUT/input.in $MY_INPUT/meshq.e  $MY_INPUT/disp.out

fwd