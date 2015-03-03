
# ==========================
# BEDTools Makefile
# (c) 2009 Aaron Quinlan
# ==========================

SHELL := /bin/bash -e

VERSION_FILE=./src/utils/version/version_git.h
RELEASED_VERSION_FILE=./src/utils/version/version_release.txt


# define our object and binary directories
export OBJ_DIR	= obj
export BIN_DIR	= bin
export SRC_DIR	= src
export UTIL_DIR	= src/utils
export SCRIPTS_DIR = scripts
export CXX		= g++
##ifeq ($(DEBUG),1)
export CXXFLAGS = -Wall -O0 -g -fno-inline -fkeep-inline-functions -D_FILE_OFFSET_BITS=64 -fPIC -DDEBUG -D_DEBUG
##else
##export CXXFLAGS = -Wall -O2 -D_FILE_OFFSET_BITS=64 -fPIC
#endif
export LIBS		= -lz
export BT_ROOT  = src/utils/BamTools/
export MKFILE_DIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

SUBDIRS = $(SRC_DIR)/lumpy

UTIL_SUBDIRS =	$(SRC_DIR)/utils/bedFile \
				$(SRC_DIR)/utils/version \
				$(SRC_DIR)/utils/gzstream \
				$(SRC_DIR)/utils/fileType \
				$(SRC_DIR)/utils/bedFilePE \
				$(SRC_DIR)/utils/genomeFile \
				$(SRC_DIR)/utils/BamTools \
				$(SRC_DIR)/utils/BamTools-Ancillary \
				$(SRC_DIR)/utils/sequenceUtilities \
				$(SRC_DIR)/utils/sqlite3


all:	lumpyexpress
	[ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR)
	[ -d $(BIN_DIR) ] || mkdir -p $(BIN_DIR)

	@echo "Building lumpy:"
	@echo "========================================================="

	@$(MAKE) --no-print-directory --directory=$(BT_ROOT) api

	@for dir in $(UTIL_SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $$dir; \
		echo ""; \
	done

	@for dir in $(SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $$dir; \
		echo ""; \
	done

lumpyexpress:
	[ -d $(BIN_DIR) ] || mkdir -p $(BIN_DIR)
	cp $(SCRIPTS_DIR)/lumpyexpress $(BIN_DIR)/lumpyexpress

	> $(BIN_DIR)/lumpyexpress.config
	@echo "LUMPY_HOME=$(MKFILE_DIR)" >> $(BIN_DIR)/lumpyexpress.config
	@echo "" >> $(BIN_DIR)/lumpyexpress.config
	@echo "LUMPY=$(MKFILE_DIR)/$(BIN_DIR)/lumpy" >> $(BIN_DIR)/lumpyexpress.config
	@echo "SAMBLASTER=`which samblaster`" >> $(BIN_DIR)/lumpyexpress.config
	@echo "SAMBAMBA=`which sambamba`" >> $(BIN_DIR)/lumpyexpress.config
	@echo "SAMTOOLS=`which samtools`" >> $(BIN_DIR)/lumpyexpress.config
	@echo "PYTHON=`which python`" >> $(BIN_DIR)/lumpyexpress.config
	@echo "" >> $(BIN_DIR)/lumpyexpress.config
	@echo "PAIREND_DISTRO=$(MKFILE_DIR)/$(SCRIPTS_DIR)/pairend_distro.py" >> $(BIN_DIR)/lumpyexpress.config
	@echo "BAMGROUPREADS=$(MKFILE_DIR)/$(SCRIPTS_DIR)/bamkit/bamgroupreads.py" >> $(BIN_DIR)/lumpyexpress.config
	@echo "BAMFILTERRG=$(MKFILE_DIR)/$(SCRIPTS_DIR)/bamkit/bamfilterrg.py" >> $(BIN_DIR)/lumpyexpress.config
	@echo "BAMLIBS=$(MKFILE_DIR)/$(SCRIPTS_DIR)/bamkit/bamlibs.py" >> $(BIN_DIR)/lumpyexpress.config

.PHONY: all

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*
	@rm -Rf $(BT_ROOT)/lib
	@rm -f $(BT_ROOT)/src/api/*.o
	@rm -f $(BT_ROOT)/src/api/internal/*.o
	@rm -Rf $(BT_ROOT)/include

.PHONY: clean
