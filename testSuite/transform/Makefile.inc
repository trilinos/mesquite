TESTTRANSFORMDIR = $(testdir)/transform

TESTTRANSFORMSRC = $(TESTTRANSFORMDIR)/main.cpp

TESTTRANSFORMEXEC = $(TESTTRANSFORMDIR)/main

all_tests: $(TESTTRANSFORMEXEC)
PREFIX =

$(TESTTRANSFORMEXEC): $(TESTTRANSFORMSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) \
	-o $@ ${TESTTRANSFORMSRC} -L$(locallibdir) -lmesquite \
        $(CONFIG_LDFLAGS) -lm
