QUALMETRICDIR = $(srcdir)/QualityMetric

QUALMETRICSRC = $(QUALMETRICDIR)/MultiplyQualityMetric.cpp \
                $(QUALMETRICDIR)/QualityMetric.cpp

QUALMETRICHEADERS = \
	$(patsubst $(QUALMETRICDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QUALMETRICDIR)/*.hpp)) \
	$(patsubst $(QUALMETRICDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QUALMETRICDIR)/*.hh))

QUALMETRICOBJ  = $(patsubst $(QUALMETRICDIR)/%.cpp, obj/%.o, $(QUALMETRICSRC)) 

ALL_MSQ_OBJ += $(QUALMETRICOBJ)
