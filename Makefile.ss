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
