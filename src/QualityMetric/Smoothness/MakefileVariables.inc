QMSMOOTHNESSDIR = $(srcdir)/QualityMetric/Smoothness

QMSMOOTHNESSSRC = $(QMSMOOTHNESSDIR)/LaplacianQualityMetric.cpp

QMSMOOTHNESSHEADERS = \
	$(patsubst $(QMSMOOTHNESSDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QMSMOOTHNESSDIR)/*.hpp)) \
	$(patsubst $(QMSMOOTHNESSDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QMSMOOTHNESSDIR)/*.hh))

QMSMOOTHNESSOBJ  = $(patsubst $(QMSMOOTHNESSDIR)/%.cpp, obj/%.o, $(QMSMOOTHNESSSRC)) 
