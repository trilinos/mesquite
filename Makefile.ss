## Default definitions for tools
CXX = CC
ARCHIVER = $(CXX) -xar -o
MAKEDEPEND = /usr/local/eng_sci/cubit/bin/makedepend -Y
LINKER = ${CXX}

## Default settings for flags
DEBUG_FLAG = -g
STD_INCLUDE_FLAG = -DUSE_STD_INCLUDES -DUSE_C_PREFIX_HEADERS
AOMD_FLAG = -D_MY_HASH_TABLE_
CXXFLAGS = ${DEBUG_FLAG} ${INCLUDE} ${STD_INCLUDE_FLAG} ${AOMD_FLAG}
SYSTEM_INCLUDE = -I/usr/include -I/opt/SUNWspro/SC5.0/include/CC -I/opt/SUNWspro/SC5.0/include/CC/std
DEPEND_FLAGS = ${SYSTEM_INCLUDE} ${DEBUG_FLAG} ${INCLUDE} ${STD_INCLUDE_FLAG}
LDFLAGS = ${DEBUG_FLAG}

## The path to the root mesquite directory.  Needs to be changed
## if you are including Makefile.customize in something other
## than the main Makefile.
MSQ_BASE_DIR = .

## Use AOMD by default.  Change this macro to link in
## something else.
TSTT_LINK = ${AOMD_TSTT_LINK}

## How to link in AOMD
AOMD_LIB_DIR = ${MSQ_BASE_DIR}/external/AOMD/lib
AOMD_TSTT_LINK = -L${AOMD_LIB_DIR} -lAOMD -lm

## How to link in MDB
MDB_LIB_DIR = ${MSQ_BASE_DIR}/external/MDB/lib
EXODUS_LIB_DIR = ${MSQ_BASE_DIR}/external/exodus/lib
MDB_TSTT_LINK = -L${MDB_LIB_DIR} -lMDB -L$(EXODUS_LIB_DIR) -lexoIIv2c -lnetcdf
