OBJ_DIR = ../../../../obj/
# -------------------
# define our includes
# -------------------
INCLUDES = -I ../htslib/htslib/

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= split.cpp split.h ReadGroup.cpp ReadGroup.h  Xam.cpp Xam.h XamHeader.cpp XamHeader.h XamReader.cpp XamReader.h XamWriter.cpp XamWriter.h XamTag.h
OBJECTS= Xam.o XamHeader.o ReadGroup.o XamReader.o XamWriter.o split.o
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

LIBRARY=libxam.a

all: $(BUILT_OBJECTS)

.PHONY: all

# $(LIBRARY): $(BUILT_OBJECTS)
# 	@echo " * linking $(LIBRARY)"
# 	ar cr $@ $^
# 	-ranlib $@

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) $(CPPFLAGS)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/*
