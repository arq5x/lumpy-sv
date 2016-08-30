OBJ_DIR = obj/
BIN_DIR = bin/
SRC_DIR = src/
CXX = g++ 
# -------------------
# define our includes
# -------------------
INCLUDES = -I htslib/htslib/
SUBDIRS = htslib/ $(SRC_DIR) 

export CXXFLAGS = -Wall -O2 -Wc++11-extensions -D_FILE_OFFSET_BITS=64 -fPIC $(INCLUDES)
BUILT_OBJECTS = $(OBJ_DIR)/*.o

all:  $(OBJ_DIR) $(BIN_DIR) $(SUBDIRS)
	@echo "- Building testing binary."
	#@$(CXX) $(CXXFLAGS) -c test.cpp -o obj/test.o $(INCLUDES)
	#@$(CXX) $(CXXFLAGS) -o $(BIN_DIR)/test $(BUILT_OBJECTS) -Ihtslib -Lhtslib htslib/libhts.a -lz 	
	@echo "done."

clean:
	rm obj/*.o htslib/*.o

.PHONY: $(SUBDIRS)
$(SUBDIRS): $(OBJ_DIR) $(BIN_DIR)
	@echo "- Building in $@"
	@$(MAKE) --no-print-directory --directory=$@
	cp htslib/*.o ../../../obj/
	cp htslib/cram/*.o ../../../obj/


