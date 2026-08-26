DEBUG ?= no

uname := $(shell uname -s)
LINKFLAGS =
ifeq ($(uname),Linux)
    LINKFLAGS = -static
endif

CC = g++
FLAGS = -O3 -std=c++11

ifeq ($(strip $(DEBUG)), yes)
	FLAGS = -O0 -std=c++11 -gdwarf-2 #-g
endif

INCL = include
LIB = lib
LIBOBJS = $(LIB)/libmaster.a

SOURCE = CenterResidue Common CreateOptions Match Parse ParseOptions ProteinStruct QueryStruct \
	 Redundancy Search SearchOptions SearchResults TargetStruct msttypes
PROGRAMS = master createPDS parsePDS
EXAMPLES = example_search

OBJECTS = $(patsubst %, objs/%.o, $(SOURCE))
BINARIES = $(patsubst %, bin/%, $(PROGRAMS))
EXAMPLEBINS = $(patsubst %, bin/%, $(EXAMPLES))

all: $(BINARIES) $(EXAMPLEBINS)
lib: $(LIBOBJS)
programs: $(BINARIES)
objects: $(OBJECTS)
examples: $(EXAMPLEBINS)

$(LIBOBJS): $(OBJECTS)
	mkdir -p lib
	ar rs $(LIBOBJS) $(OBJECTS)
$(OBJECTS): objs/%.o : src/%.cpp include/%.h
	$(CC) $(FLAGS) -I$(INCL) -c $< -o $@
$(BINARIES): bin/% : programs/%.cpp $(LIBOBJS)
	$(CC) $(FLAGS) -I$(INCL) $< -L$(LIB) -lmaster $(LINKFLAGS) -o $@
$(EXAMPLEBINS): bin/% : examples/%.cpp $(LIBOBJS)
	$(CC) $(FLAGS) -I$(INCL) $< -L$(LIB) -lmaster $(LINKFLAGS) -o $@

.PHONY: clean
clean:
	-rm -f $(LIBOBJS) $(OBJECTS) $(BINARIES) $(EXAMPLEBINS)
