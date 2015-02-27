UTILITIES_DIR = ../utils/
OBJ_DIR = ../../obj/
BIN_DIR = ../../bin/

# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/bedFile/ \
           -I$(UTILITIES_DIR)/bedFilePE/ \
           -I$(UTILITIES_DIR)/genomeFile/ \
           -I$(UTILITIES_DIR)/version/ \
           -I$(UTILITIES_DIR)/gzstream/ \
           -I$(UTILITIES_DIR)/lineFileUtilities/ \
           -I$(UTILITIES_DIR)/fileType/ \
           -I$(UTILITIES_DIR)/BamTools/include \
		   -I$(UTILITIES_DIR)/BamTools-Ancillary \
           -I$(UTILITIES_DIR)/sequenceUtilities \
           -I$(UTILITIES_DIR)/UCSCBins \
           -I$(UTILITIES_DIR)/sqlite3 

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= lumpy.cpp \
		 SV_Pair.cpp \
		 SV_Evidence.cpp \
		 SV_BreakPoint.cpp \
		 SV_Bedpe.cpp \
		 SV_SplitRead.cpp \
		 SV_SplitReadReader.cpp \
		 SV_PairReader.cpp \
		 SV_BedpeReader.cpp \
		 SV_BamReader.cpp \
		 SV_InterChromBamReader.cpp \
		 SV_EvidenceReader.cpp \
		 SV_Tools.cpp \
		 SV_VcfVariant.cpp \
		 log_space.cpp
OBJECTS= $(SOURCES:.cpp=.o)
_EXT_OBJECTS=BamAncillary.o \
			 sequenceUtils.o \
			 gzstream.o \
			 bedFilePE.o \
			 bedFile.o \
			 genomeFile.o \
			 fileType.o \
			 sqlite3.o 
EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
PROGRAM=lumpy


all: $(PROGRAM)

$(PROGRAM): $(BUILT_OBJECTS) $(EXT_OBJECTS)
	@echo "  * linking $(PROGRAM)"
	@$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/$@ $^ \
		-L$(UTILITIES_DIR)/BamTools/lib/ -lbamtools \
		$(LIBS) \
		-lz -lm -ldl -lpthread

$(OBJ_DIR)/%.o: %.cpp
	@echo "  * compiling" $<
	@$(CXX) -g  -c $(CXXFLAGS) $(LDFLAGS) $(INCLUDES) $< -o $@

$(EXT_OBJECTS):
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/BamTools-Ancillary/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/BamTools/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/bedFilePE/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/sequenceUtilities/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/lineFileUtilities/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/gzstream/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/fileType/ 
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/sqlite3/ 

clean:
	@echo "Cleaning up."
	rm -f $(BUILT_OBJECTS)

.PHONY: clean
