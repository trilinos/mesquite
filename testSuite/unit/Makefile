# Change this in Makefile.customize so that it points to your
# installation of cppunit
#CPPUNIT_DIR = ../../../cppunit


# The rest of these items probably won't need to be messed with
#CPPUNIT_INCLUDE = -I${CPPUNIT_DIR}/include
MSQ_BASE_DIR = ../..
TSTT_INCLUDE = -I$(MSQ_BASE_DIR)/TSTT-interface
MESQUITE_INCLUDE = ${TSTT_INCLUDE} -I$(MSQ_BASE_DIR)/include -I$(MSQ_BASE_DIR)/includeLinks
TEST_INCLUDE = ${CPPUNIT_INCLUDE} ${MESQUITE_INCLUDE}

MESQUITE_LINK = -L$(MSQ_BASE_DIR)/lib -lmesquite ${TSTT_LINK}
MESQUITE_AOMD_LINK = -L$(MSQ_BASE_DIR)/lib -lmesquite ${AOMD_TSTT_LINK}
MESQUITE_MDB_LINK = -L$(MSQ_BASE_DIR)/lib -lmesquite ${MDB_TSTT_LINK}

TEST_LINK = ${CPPUNIT_LINK} ${MESQUITE_LINK} -lm
TEST_AOMD_LINK = ${CPPUNIT_LINK} ${MESQUITE_AOMD_LINK} -lm
TEST_MDB_LINK = ${CPPUNIT_LINK} ${MESQUITE_MDB_LINK} -lm
PREFIX = @
PREFIX = 

default: msq_test

include Makefile.inc
include $(MSQ_BASE_DIR)/Makefile.customize
MSQ_BASE_DIR = ../..

msq_test : $(TEST_OBJS) $(MSQ_BASE_DIR)/lib/libmesquite.a
	@echo "Linking $@"
	$(PREFIX)$(LINKER) $(LDFLAGS) \
        $(TEST_OBJS) \
        -o $@ \
        -L$(MSQ_BASE_DIR)/lib -lmesquite \
        -L$(CPPUNIT_DIR)/lib -lcppunit -ldl

msq_test_aomd : ${TEST_OBJS}
	@echo "Linking $@"
	${PREFIX}${LINKER} ${LDFLAGS} ${TEST_OBJS} ${TEST_AOMD_LINK} -o $@

msq_test_mdb : ${TEST_OBJS}
	@echo "Linking $@"
	${PREFIX}${LINKER} ${LDFLAGS} ${TEST_OBJS} ${TEST_MDB_LINK} -o $@

clean:
	@rm -f *.o msq_test
	@rm -f *.o msq_test_mdb
	@rm -f *.o msq_test_aomd
	@rm -rf SunWS_cache
	@rm -rf cxx_repository

.SUFFIXES : .cpp .o
