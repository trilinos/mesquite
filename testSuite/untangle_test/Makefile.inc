UNTDIR = $(testdir)/untangle_test

UNTSRC = $(UNTDIR)/main.cpp

UNTEXEC = $(UNTDIR)/main

all_tests: $(UNTEXEC)

$(UNTEXEC): $(UNTSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) \
	-o $@ ${UNTSRC} -L$(locallibdir) -lmesquite \
        $(CONFIG_LDFLAGS) -lm

