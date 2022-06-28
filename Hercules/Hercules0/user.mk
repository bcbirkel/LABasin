SYSTEM = USCHPC
RUN_DIR = /home/scec-00/birkel/Hercules0/quake/forward/

#CFLAGS += -g -O0 -DDEBUG

CFLAGS += -O3 -msse2
CFLAGS += -DCVM_SRCPATH=\"$(SCRATCH)/cvmetrees/terav2.e\"
CFLAGS += -DCVM_DESTDIR=\"$(RUN_DIR)/\"

#LDFLAGS += -g

# Vis
#VIS_CFLAGS = -DVIS -DVISINTERVAL=10
VIS_CFLAGS =

# I/O
# prevent database replication
# definining the SCEC macro prevents CVM database replication

IO_CFLAGS = -DUSECVMDB -DSCEC -DPROCPERNODE=4000
#IO_CFLAGS -DNO_OUTPUT


