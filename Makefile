
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
export CXX		= g++
ifeq ($(DEBUG),1)
export CXXFLAGS = -Wall -O0 -g -fno-inline -fkeep-inline-functions -D_FILE_OFFSET_BITS=64 -fPIC -DDEBUG -D_DEBUG
else
export CXXFLAGS = -Wall -O2 -D_FILE_OFFSET_BITS=64 -fPIC
endif
export LIBS		= -lz
export BT_ROOT  = src/utils/BamTools/


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


all:
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


.PHONY: all

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*
	@rm -Rf $(BT_ROOT)/lib
	@rm -f $(BT_ROOT)/src/api/*.o
	@rm -f $(BT_ROOT)/src/api/internal/*.o
	@rm -Rf $(BT_ROOT)/include

.PHONY: clean
