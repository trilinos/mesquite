QMUNTANGLEDIR = $(srcdir)/QualityMetric/Untangle

QMUNTANGLESRC = 

QMUNTANGLEHEADERS = \
	$(patsubst $(QMUNTANGLEDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QMUNTANGLEDIR)/*.hpp)) \
	$(patsubst $(QMUNTANGLEDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QMUNTANGLEDIR)/*.hh))

QMUNTANGLEOBJ  = $(patsubst $(QMUNTANGLEDIR)/%.cpp, obj/%.o, $(QMUNTANGLESRC)) 
