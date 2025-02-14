# Load the configuration file
include ../config.mk
include ./Make.defs

# RMT3D requires PGMG library
# Input PGMG path in Make.defs, otherwise, default is ../../
PGMG_DIR?=../../
# Assemble the subdirectories that contain source files
source_dirs += ../F3D_Source ../Extrapolation ../IO ../Sim ../Utilities

# source files and corresponding objects
src+=$(foreach subdir,$(source_dirs),$(wildcard $(subdir)/*.cc))
objs=$(patsubst %.cc,%.o,$(src))

# header files
headers+=$(foreach subdir,$(source_dirs),$(wildcard $(subdir)/*.hh))

# Add these directory to include path
iflags = $(foreach subdir,$(source_dirs), -I./$(subdir))
# Add a few flags for GNU Scientific Library
iflags+=-I$(PGMG_DIR) `gsl-config --cflags`
# Add these directory to library path
lflags+=-L. -L$(PGMG_DIR) -lpgmg `gsl-config --libs`
# add C++-11 std to clfags

########################################
# Definition of macros in Make.defs
ifeq ($(DEBUG), TRUE)
	cflags += -DDEBUG
	debugsuffix := .debug
else
	debugsuffix :=
endif

# COMPILE A VARIABLE DENSITY VERSION
ifeq ($(VAR_DEN), TRUE)
	cflags += -DVAR_DEN
	vardensuffix := .varden
else
	vardensuffix :=
endif

# COMPILE A VERSION THAT DOES TANGENTIAL STABILITY CALCULATION
ifeq ($(TANG), TRUE)
	cflags += -DTANG
	tangsuffix := .tang
else
	tangsuffix :=
endif


execsuffix :=$(vardensuffix)$(debugsuffix)$(tangsuffix)

########################################
# required executables
execs=ctr_stats unpack tr_unpack run_sim$(execsuffix)

########################################
# make everything in pgmg and then here
all:
	$(MAKE) depend
	$(MAKE) -C $(PGMG_DIR) objects
	$(MAKE) -C ../pov_headers
	$(MAKE) executables

executables: $(execs)

-include Makefile.dep

depend: $(src)
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

clean:
	rm -rf $(execs) $(objs) librmt3d.a *.dSYM

realclean:
	rm -rf $(execs) $(objs) librmt3d.a *.dSYM
	$(MAKE) -C $(PGMG_DIR)/ clean


# object rules
%.o: %.cc
	$(mpicxx) $(cflags) $(cpp11) $(iflags) -c $< -o $@

librmt3d.a: $(objs)
	rm -f librmt3d.a
	ar rs librmt3d.a $^

run_sim$(execsuffix): run_sim.cc librmt3d.a $(headers)
	$(mpicxx) $(cflags) $(iflags) -o $@ $^ $(lflags) -lrmt3d $(lp_lflags)

ctr_stats: ctr_stats.cc librmt3d.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) -lrmt3d $(lp_lflags) -g -O0

unpack: unpack.cc librmt3d.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) -lrmt3d $(lp_lflags)

tr_unpack: tr_unpack.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

.PHONY: clean all executables depend realclean
