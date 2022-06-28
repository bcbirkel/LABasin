#!/bin/bash

# RUN HERCURLES SIMULATIONS IN HPC-USC
# GIVE EQK NAME, SIMULATION NUMBER, FORWARD|ADJOINT, SOURCE NUMBER, SIM TIME & NUMBER OF PROCS
EQKN=B20140329
SIMN=0
FOAD=Forward
SRCN=Source1
TIME=10:00:00
NTSK=16

# HERCULES IS CHOSEN BASED ON SIM NUMBER, IF REPEATED MAYBE AND ERROR
let "HERN=$SIMN%10"
echo "RUNNING $EQKN SIM$SIMN - $FOAD - SRC$SRCN IN HERCULES $HERN"

# COPY SLURM FILE TO SIMULATION NUMBER
cp Forward.sl Forward${SIMN}.sl

# SUBMIT SLURM JOB (REVIEW PARAMETERS IN FILE.sl)
# TAKEN OUT --partition=scec \
sbatch --export=E=$EQKN,F=$FOAD,T=$NTSK,H=$HERN,S=$SIMN,R=$SRCN \
	--ntasks=$NTSK --time=$TIME \
	--output=${EQKN}Sim${SIMN}.o \
	--error=${EQKN}Sim${SIMN}.e \
	--constraint=IB \
	--partition=scec \
	--job-name=Sim${SIMN} Forward${SIMN}.sl 
	
# STAY CONNECTED AND CHECKING JOB STATUS
SQUEUE=$(squeue -u birkel | grep "Sim")
while [ -n "$SQUEUE" ]; do
        squeue -u birkel
        echo "RUNNING SIMULATION $SIMN"
       	sleep 2m
       	SQUEUE=$(squeue -u birkel | grep "Sim")
done
	

