UNTDIR = $(testdir)/untangle_test

UNTSRC = $(UNTDIR)/main.cpp

UNTEXEC = $(UNTDIR)/main

all_tests: $(UNTEXEC)

$(UNTEXEC): $(UNTSRC) $(locallibdir)/mesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} -o $@ ${UNTSRC} $(locallibdir)/mesquite.a -Lexternal/AOMD/lib/ -lAOMD

