#  MESQUITE -- The Mesh Quality Improvement Toolkit
#
#  Copyright 2004 Lawrence Livermore National Laboratory.  Under 
#  the terms of Contract B545069 with the University of Wisconsin -- 
#  Madison, Lawrence Livermore National Laboratory retains certain
#  rights in this software.
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License 
#  (lgpl.txt) along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#  kraftche@cae.wisc.edu    

# This directory
TSTT_DIR = tstt

# Subdirectories to create
TSTT_REPO_DIR = $(TSTT_DIR)/repo
TSTT_SHAR_DIR = $(TSTT_DIR)/share
TSTT_MESH_DIR = $(TSTT_DIR)/mesh
TSTT_GEOM_DIR = $(TSTT_DIR)/geom
TSTT_ASSO_DIR = $(TSTT_DIR)/assoc

# Input files for Babel
TSTT_SHAR_SIDL = $(TSTT_DIR)/TSTT.sidl
TSTT_MESH_SIDL = $(TSTT_DIR)/TSTTM.sidl
TSTT_GEOM_SIDL = $(TSTT_DIR)/TSTTG.sidl
TSTT_ASSO_SIDL = $(TSTT_DIR)/TSTTR.sidl


# Include defines for compiling TSTT-related code
TSTT_SHAR_INC = -I$(TSTT_SHAR_DIR)
TSTT_GEOM_INC = -I$(TSTT_GEOM_DIR) -I$(TSTT_ASSO_DIR)
TSTT_MESH_INC = -I$(TSTT_MESH_DIR)
ifeq ($(MSQ_TSTT_MESH),yes)
  TSTT_INC += $(TSTT_MESH_INC) $(TSTT_SHAR_INC)
endif
ifeq ($(MSQ_TSTT_GEOM),yes)
  TSTT_INC += $(TSTT_GEOM_INC) $(TSTT_SHAR_INC)
endif

# Utility library names for linking test code
TSTT_SHAR_CLIENT_LIBNAME = msqtsttclient
TSTT_GEOM_CLIENT_LIBNAME = msqtsttgeom
TSTT_MESH_CLIENT_LIBNAME = msqtsttmesh

# Link defines for link TSTT-related tests
TSTT_SHAR_LNK = -l$(TSTT_SHAR_CLIENT_LIBNAME)
TSTT_GEOM_LNK = -l$(TSTT_GEOM_CLIENT_LIBNAME)
TSTT_MESH_LNK = -l$(TSTT_MESH_CLIENT_LIBNAME)

# Library names for TSTT-related test dependencies
TSTT_SHAR_CLIENT_LIB = $(locallibdir)/lib$(TSTT_SHAR_CLIENT_LIBNAME).a
TSTT_GEOM_CLIENT_LIB = $(locallibdir)/lib$(TSTT_GEOM_CLIENT_LIBNAME).a
TSTT_MESH_CLIENT_LIB = $(locallibdir)/lib$(TSTT_MESH_CLIENT_LIBNAME).a


# Filename to use as a time stamp
# Don't change this.  GNU Make has an implicit dependency
# on the specified files if they're included.  This is
# what causes the headers and stubs to be generated from
# the .sidl files.
TSTT_TARGET_FILE  = babel.make

# Construct short-hand target names from above variables
TSTT_REPO_TARGET = $(TSTT_REPO_DIR)/.timestamp
TSTT_SHAR_TARGET = $(TSTT_SHAR_DIR)/$(TSTT_TARGET_FILE)
TSTT_MESH_TARGET = $(TSTT_MESH_DIR)/$(TSTT_TARGET_FILE)
TSTT_GEOM_TARGET = $(TSTT_GEOM_DIR)/$(TSTT_TARGET_FILE)
TSTT_ASSO_TARGET = $(TSTT_ASSO_DIR)/$(TSTT_TARGET_FILE)


# Rules for generating headers and source stubs from .sidl files

# Build XML repository for shared TSTT definitions
# Used by other TSTT interface definitions that depend 
# on the shared interface.
$(TSTT_REPO_TARGET): $(TSTT_SHAR_SIDL)
	@echo "Generating shared TSTT interface XML repository"
	$(PREFIX) rm -f $(TSTT_REPO_DIR)/* 
	$(PREFIX) $(BABEL) -tXML -o$(TSTT_REPO_DIR) $<
	$(PREFIX) touch $(TSTT_REPO_TARGET)

# Build headers for shared TSTT defintions
$(TSTT_SHAR_TARGET): $(TSTT_SHAR_SIDL) $(TSTT_REPO_TARGET)
	@echo "Generating shared TSTT interface headers and stubs"
	$(PREFIX) rm -f $(TSTT_SHAR_DIR)/*
	$(PREFIX) $(BABEL) -cC++ -o$(TSTT_SHAR_DIR) $<

# Build headers for TSTT mesh interface
$(TSTT_MESH_TARGET): $(TSTT_MESH_SIDL) $(TSTT_SHAR_TARGET)
	@echo "Generating TSTT mesh interface headers and stubs"
	$(PREFIX) rm -f $(TSTT_MESH_DIR)/*
	$(PREFIX) $(BABEL) -cC++ -E -R$(TSTT_REPO_DIR) -o$(TSTT_MESH_DIR) $<

# Build headers for TSTT geometry interface
$(TSTT_GEOM_TARGET): $(TSTT_GEOM_SIDL) $(TSTT_SHAR_TARGET)
	@echo "Generating TSTT geometry interface headers and stubs"
	$(PREFIX) rm -f $(TSTT_GEOM_DIR)/*
	$(PREFIX) $(BABEL) -cC++ -E -R$(TSTT_REPO_DIR) -o$(TSTT_GEOM_DIR) $<

# Build headers for TSTT mesh-geometry association interface
$(TSTT_ASSO_TARGET): $(TSTT_ASSO_SIDL) $(TSTT_SHAR_TARGET)
	@echo "Generating TSTT classification interface headers and stubs"
	$(PREFIX) rm -f $(TSTT_ASSO_DIR)/*
	$(PREFIX) $(BABEL) -cC++ -E -R$(TSTT_REPO_DIR) -o$(TSTT_ASSO_DIR) $<



# Rules to build temporary libaries containing TSTT stub code for
# linking TSTT tests

ifeq ($(MSQ_TSTT_MESH),yes)
  include $(TSTT_MESH_TARGET)
  TSTT_MESH_STUBS := $(STUBSRCS:%=$(TSTT_MESH_DIR)/%)
  MSQ_TSTT_COMMON = 1
  TSTT_STUB_LIBS += $(TSTT_MESH_CLIENT_LIB)
endif

ifeq ($(MSQ_TSTT_GEOM),yes)
  include $(TSTT_GEOM_TARGET)
  TSTT_GEOM_STUBS := $(STUBSRCS:%=$(TSTT_GEOM_DIR)/%)
  include $(TSTT_ASSO_TARGET)
  TSTT_ASSO_STUBS := $(STUBSRCS:%=$(TSTT_ASSO_DIR)/%)
  MSQ_TSTT_COMMON = 1
  TSTT_STUB_LIBS += $(TSTT_GEOM_CLIENT_LIB)
endif

ifdef MSQ_TSTT_COMMON
  include $(TSTT_SHAR_TARGET)
  TSTT_SHAR_STUBS := $(STUBSRCS:%=$(TSTT_SHAR_DIR)/%)
  TSTT_STUB_LIBS += $(TSTT_SHAR_CLIENT_LIB)
endif


all_libs: $(TSTT_STUB_LIBS)

$(TSTT_GEOM_CLIENT_LIB): $(TSTT_GEOM_STUBS:.cc=.o)
	$(PREFIX) $(ARCHIVER) $@ $^

$(TSTT_MESH_CLIENT_LIB): $(TSTT_MESH_STUBS:.cc=.o)
	$(PREFIX) $(ARCHIVER) $@ $^

$(TSTT_SHAR_CLIENT_LIB): $(TSTT_SHAR_STUBS:.cc=.o)
	$(PREFIX) $(ARCHIVER) $@ $^

# Compiling rules
# Mesquite uses .cpp files, so we should be able to define
# a rule here for .cc files w/out messing up the rest of
# the Mesquite build
.cc.o:
	$(PREFIX) $(CXX) -c $< -o $@ $(CXXFLAGS) \
	$(SIDL_INC) $(TSTT_SHAR_INC) $(TSTT_GEOM_INC) $(TSTT_MESH_INC)

# Misc rules

.PHONEY: tstt_clean tstt_veryclean
clean: tstt_clean
veryclean: tstt_veryclean
tstt_clean: 
	$(PREFIX) rm -rf $(TSTT_REPO_DIR)
	$(PREFIX) rm -rf $(TSTT_SHAR_DIR)
	$(PREFIX) rm -rf $(TSTT_MESH_DIR)
	$(PREFIX) rm -rf $(TSTT_GEOM_DIR)
	$(PREFIX) rm -rf $(TSTT_ASSO_DIR)
tstt_veryclean: tstt_clean
	$(PREFIX) rm -f $(TSTT_SHAR_CLIENT_LIB)
	$(PREFIX) rm -f $(TSTT_MESH_CLIENT_LIB)
	$(PREFIX) rm -f $(TSTT_GEOM_CLIENT_LIB)
