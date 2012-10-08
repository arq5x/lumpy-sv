UTILITIES_DIR = ../utils/
OBJ_DIR = ../../obj/
BIN_DIR = ../../bin/

GSL_INCLUDE=-I/sw/include/ \
			-I/sw/include/gsl
GSL_LINK=-L/sw/lib/

# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/bedFile/ \
           -I$(UTILITIES_DIR)/bedFilePE/ \
           -I$(UTILITIES_DIR)/version/ \
           -I$(UTILITIES_DIR)/gzstream/ \
           -I$(UTILITIES_DIR)/lineFileUtilities/ \
           -I$(UTILITIES_DIR)/fileType/ \
           -I$(UTILITIES_DIR)/BamTools/include \
		   -I$(UTILITIES_DIR)/BamTools-Ancillary \
           -I$(UTILITIES_DIR)/sequenceUtilities/ \
           -I$(UTILITIES_DIR)/UCSCBins \
           $(GSL_INCLUDE)

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
		 SV_EvidenceReader.cpp \
		 SV_Tools.cpp \
		 log_space.cpp
OBJECTS= $(SOURCES:.cpp=.o)
_EXT_OBJECTS=BamAncillary.o \
			 sequenceUtils.o \
			 gzstream.o \
			 bedFilePE.o \
			 fileType.o 
EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
PROGRAM=lumpy


all: $(PROGRAM)

$(PROGRAM): $(BUILT_OBJECTS) $(EXT_OBJECTS)
	@echo "  * linking $(PROGRAM)"
	@$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/$@ $^ \
		$(GSL_LINK) -L$(UTILITIES_DIR)/BamTools/lib/ -lbamtools \
		$(LIBS) \
		-lz -lgsl  -lm

$(OBJ_DIR)/%.o: %.cpp
	@echo "  * compiling" $<
	@$(CXX) -g  -c $(CXXFLAGS) $(LDFLAGS) $(INCLUDES) $< -o $@

$(EXT_OBJECTS):
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/BamTools-Ancillary/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/bedFilePE/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/sequenceUtilities/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/lineFileUtilities/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/gzstream/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/fileType/ 
	
clean:
	@echo "Cleaning up."
	rm -f $(BUILT_OBJECTS) $(EXT_OBJECTS)

.PHONY: clean
