TESTG843DIR = $(testdir)/Guides843

TESTG843SRC = $(TESTG843DIR)/main.cpp

TESTG843EXEC = $(TESTG843DIR)/main

all_tests: $(TESTG843EXEC)


$(TESTG843EXEC): $(TESTG843SRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} -o $@ ${TESTG843SRC} -L$(locallibdir) -lmesquite ${TSTT_LINK} -lm
