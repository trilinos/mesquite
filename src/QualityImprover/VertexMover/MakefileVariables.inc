VERTEXMOVERDIR = $(srcdir)/QualityImprover/VertexMover

VERTEXMOVERSRC = $(VERTEXMOVERDIR)/VertexMover.cpp 


VERTEXMOVERHEADERS = \
	$(patsubst $(VERTEXMOVERDIR)/%,$(localincludedir)/%,\
	$(wildcard $(VERTEXMOVERDIR)/*.hpp)) \
	$(patsubst $(VERTEXMOVERDIR)/%,$(localincludedir)/%,\
	$(wildcard $(VERTEXMOVERDIR)/*.hh))

VERTEXMOVEROBJ  = $(patsubst $(VERTEXMOVERDIR)/%.cpp, obj/%.o, $(VERTEXMOVERSRC)) 

ALL_MSQ_OBJ += $(VERTEXMOVEROBJ)
