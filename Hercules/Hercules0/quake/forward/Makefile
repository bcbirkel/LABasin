# -*- Makefile -*-
# 
#
ifndef WORKDIR
	WORKDIR = $(PWD)/../..
endif

include $(WORKDIR)/systemdef.mk
include $(WORKDIR)/common.mk

#CFLAGS += -DDEBUG_OVERFLOW
#CFLAGS += -DAUTO_DELTA_T
#CFLAGS += -DVALIDATION

CFLAGS += -I$(ETREE_DIR) -I$(CVM_DIR) -I$(OCTOR_DIR) -I$(VIS_DIR) -I.
CFLAGS += $(SOLVE_CFLAGS) $(IO_CFLAGS) $(VIS_CFLAGS)
LOADLIBES = $(CVM_DIR)/cvm.o $(OCTOR_LIB) $(ETREE_LIB) 

ifeq ($(SYSTEM), BGW)
	LOADLIBES += /bgl/local/lib/libmpihpm_c.a /bgl/BlueLight/ppcfloor/bglsys/lib/libbgl_perfctr.rts.a
endif

#
# Object modules 
#
OBJECTS = qmesh.o psolve.o q4.o single_query.o nrutila.o quakesource.o geometrics.o geodataqinit.o geodataqquery.o  geominterp.o Stdio.o
.PHONY: all tags clean cleanall move_cxx_repository

all: move_cxx_repository $(TARGET)

TARGET = qmesh psolve single_query q4showmeta

move_cxx_repository:
ifeq ($(SYSTEM), LEMIEUX)
	cp -rf $(VIS_DIR)/cxx_repository .
endif

mesh_query: mesh_query.o query_mesh_property.o Stdio.o

qmesh: qmesh.o

psolve: psolve.o nrutila.o quakesource.o geometrics.o output.o util.o  geometrics.o geodataqinit.o geodataqquery.o geominterp.o Stdio.o
	$(LD) $^ \
	$(LOADLIBES) \
	$(LDFLAGS) \
	-o $@

q4showmeta: q4showmeta.o util.o

psolve_output: psolve_output.o

single_query: single_query.o q4.o

output.o: output.c output.h psolve.h
psolve.o: psolve.c Stdio.c psolve.h output.h quakesource.h
quakesource.o: quakesource.c quakesource.h psolve.h geodataqinit.h
geodataqinit.o: geodataqinit.c Stdio.c Stdio.h geodataqinit.h
Stdio.o: Stdio.c Stdio.h
geodataqquery.o: geodataqquery.c Stdio.c geodataq.h geominterp.h Stdio.h
geominterp.o: geominterp.c geominterp.h geodataq.h


all: $(TARGET)

clean:
	rm -rf $(OBJECTS)  core *.o *~ ./cxx_repository

cleanall:
	rm -rf $(OBJECTS)  $(TARGET) core *~ ./cxx_repository

TAGS:
	etags *.c *.h

tags:	TAGS

# $Id: Makefile,v 1.8 2007/07/12 23:55:22 jclopez Exp $
