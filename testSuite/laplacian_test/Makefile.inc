LAPLDIR = $(testdir)/laplacian_test

LAPLSRC = $(LAPLDIR)/main.cpp

LAPLEXEC = $(LAPLDIR)/main

all_tests: $(LAPLEXEC)

$(LAPLEXEC): $(LAPLSRC) $(locallibdir)/mesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} -o $@ ${LAPLSRC} $(locallibdir)/mesquite.a -Lexternal/AOMD/lib/ -R../../external/AOMD/lib/ -lAOMD
