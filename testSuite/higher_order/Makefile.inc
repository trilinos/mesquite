HIGHERORDERDIR = $(testdir)/higher_order

HIGHERORDERSRC = $(HIGHERORDERDIR)/main.cpp

HIGHERORDEREXEC = $(HIGHERORDERDIR)/main

all_tests: $(HIGHERORDEREXEC)

$(HIGHERORDEREXEC): $(HIGHERORDERSRC) $(locallibdir)/libmesquite.a
	@echo "Linking Mesquite driver application ($@)..."
	${LINKER} $(CXXFLAGS) ${LDFLAGS} $(CONFIG_CFLAGS) \
	-o $@ ${HIGHERORDERSRC} -L$(locallibdir) -lmesquite \
        $(CONFIG_LDFLAGS) -lm
