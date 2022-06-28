# -*- Makefile -*-
#
# Description:	Quake's Hercules toolchain top directory Makefile
#
-include user.mk

WORKDIR = $(CURDIR)

#ifndef TARGET_DIR
	TARGET_DIR = /lustre/ttu/r1/1024.image
#endif

include systemdef.mk
include common.mk

.PHONY: all clean anim

all:
	$(MAKE) -C $(ETREE_DIR)   SYSTEM=$(SYSTEM) WORKDIR=$(WORKDIR)
	$(MAKE) -C $(OCTOR_DIR)   SYSTEM=$(SYSTEM) WORKDIR=$(WORKDIR)
	$(MAKE) -C $(CVM_DIR)     SYSTEM=$(SYSTEM) WORKDIR=$(WORKDIR)
	$(MAKE) -C $(FORWARD_DIR) SYSTEM=$(SYSTEM) WORKDIR=$(WORKDIR)

clean:
	$(MAKE) -C $(ETREE_DIR)   WORKDIR=$(WORKDIR) clean
	$(MAKE) -C $(OCTOR_DIR)   WORKDIR=$(WORKDIR) clean
	$(MAKE) -C $(CVM_DIR)     WORKDIR=$(WORKDIR) clean
	$(MAKE) -C $(FORWARD_DIR) WORKDIR=$(WORKDIR) clean


cleanall: MY_DIRS= $(ETREE_DIR) $(OCTOR_DIR) $(CVM_DIR) $(FORWARD_DIR)
cleanall: clean
	for MY_DIR in $(MY_DIRS) ; \
	   do $(MAKE) -C $${MY_DIR} WORKDIR=$(WORKDIR) cleanall ; \
	done

# $Id: Makefile,v 1.5 2006/07/31 04:24:43 jclopez Exp $
