# -*- Makefile -*-
# Override SYSTEM and other variable definition in user.mk
# In order to switch between a debug and an optimized executable, set CFLAGS
# in your 'user.mk' file as follows:
# * for an executable with debug information:
#	CFLAGS = -g -DDEBUG -O0
# * for an optimized executable:
#	CFLAGS = -O2
# check other platform specific flags below.

-include $(WORKDIR)/user.mk

ifndef SYSTEM
	SYSTEM = USCHPC
endif


ifeq ($(SYSTEM), BGW)
	MPI_DIR=/bgl/BlueLight/ppcfloor/bglsys
	CC	= mpcc
	CXX	= mpCC
	LD	=  mpCC
	CFLAGS += -qmaxmem=-1 -qarch=440
	CFLAGS += -DBGL 
	CFLAGS += -D_FILE_OFFSET_BITS=64 
endif


ifeq ($(SYSTEM), BGL)
	MPI_DIR=/bgl/BlueLight/ppcfloor/bglsys
	CC	= mpicc
	CXX	= mpicxx
	LD	= mpicxx
	CFLAGS += -DBGL 
	CFLAGS += -D_FILE_OFFSET_BITS=64 -Wall -D_LARGEFILE_SOURCE 
endif

# Note: Add the following flags ( -fastsse -O3 ) to the CFLAGS variable
# in order to create an optimized executable.
#
# Also define USE_IOBUF=1 in order to use the io_buf library and remember
# to load the io_buf module before compiling and launching your job.

ifeq ($(SYSTEM), BIGBEN)
	CC 	= cc
	CXX 	= CC
	LD 	= CC
	CFLAGS  += -DBIGBEN -target=catamount
	LDFLAGS += -target=catamount
	ifdef IOBUF_INC
	    CPPFLAGS += -I${IOBUF_INC}
	endif
endif

# optimization flags: -fast -O
# debug flags: 	-DDEBUG -g

ifeq ($(SYSTEM), LEMIEUX)
	CC	= cc
	CXX	= cxx
	LD	= cxx
	CFLAGS  += -DPROCPERNODE=4
	CFLAGS  += -trapuv -check_bounds -DALPHA_TRU64UNIX_CC -DALIGNMENT
	LDFLAGS += -lelan -lmpi
endif


ifeq ($(SYSTEM), MANTEO)
	MPI_DIR	= /usr/local/mpich-1.2.6
	CC	= $(MPI_DIR)/bin/icc
	CXX	= $(MPI_DIR)/bin/icpc
	LD	= $(MPI_DIR)/bin/icpc
	CFLAGS += -Wall -DPROCPERNODE=64
	CFLAGS += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
	CFLAGS += -DCVM_SRCPATH=\"/usr0/tutk/cvmetrees/labase.e\"
	CFLAGS += -DCVM_DESTDIR=\"/usr1/tutk/cvmetrees\"
endif

# Add, define or override the following options in $(WORKDIR)/user.mk:
#   RUN_DIR = /path/to/your/run/time/directory
#   CFLAGS += -g
#   CFLAGS += -g -DDEBUG
#   CFLAGS += -DPROCPERNODE=2
#   CFLAGS += -DCVM_SRCPATH=\"/path/to/your/cvm/database/file\"
#   CFLAGS += -DCVM_DESTDIR=\"$(RUN_DIR)/destination/directory\"
#
ifeq ($(SYSTEM), LINUX-MPICH)
	RUN_DIR =\"/home/scec-00/birkel/Software/Hercules/quake/forward/\"
	MPI_DIR =/usr/usc/openmpi/default/
	CFLAGS += -DPROCPERNODE=2
	CFLAGS += -DCVM_SRCPATH=\"/home/scec-00/birkel/Software/Hercules/quake/forward/labasin.e\"
	CFLAGS += -DCVM_DESTDIR=\"$(RUN_DIR)/\"
	CC      = $(MPI_DIR)/bin/mpicc
	CXX     = $(MPI_DIR)/bin/mpicc
	LD      = $(MPI_DIR)/bin/mpicc
	CPPFLAGS += -DMPICH_IGNORE_CXX_SEEK
	CFLAGS += -Wall -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
endif


ifeq ($(SYSTEM), USCHPC)
        MPI_DIR = /usr/usc/openmpi/1.8.7/
        CC	= $(MPI_DIR)/bin/mpicc
	CXX	= $(MPI_DIR)/bin/mpicxx
        LD	= $(MPI_DIR)/bin/mpicxx
        CFLAGS += -DPROCPERNODE=1
        CFLAGS += -Wall -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
        CFLAGS += -DCVM_SRCPATH=\"/home/scec-00/birkel/Hercules0/quake/forward/labase.e\"
        CFLAGS += -DCVM_DESTDIR=\"/tmp\"
endif

