QUALIMPROVERDIR = $(srcdir)/QualityImprover

QUALIMPROVERSRC = 

QUALIMPROVERHEADERS = \
	$(patsubst $(QUALIMPROVERDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QUALIMPROVERDIR)/*.hpp)) \
	$(patsubst $(QUALIMPROVERDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QUALIMPROVERDIR)/*.hh))

QUALIMPROVEROBJ  = $(patsubst $(QUALIMPROVERDIR)/%.cpp, obj/%.o, $(QUALIMPROVERSRC))
