QMVOLUMEDIR = $(srcdir)/QualityMetric/Volume

QMVOLUMESRC = 

QMVOLUMEHEADERS = \
	$(patsubst $(QMVOLUMEDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QMVOLUMEDIR)/*.hpp)) \
	$(patsubst $(QMVOLUMEDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QMVOLUMEDIR)/*.hh))

QMVOLUMEOBJ  = $(patsubst $(QMVOLUMEDIR)/%.cpp, obj/%.o, $(QMVOLUMESRC))
