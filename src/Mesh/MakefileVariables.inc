MESHDIR = $(srcdir)/Mesh

MESHSRC = $(MESHDIR)/MeshSet.cpp \
	       $(MESHDIR)/MsqMeshEntity.cpp \
	       $(MESHDIR)/MsqVertex.cpp \
	       $(MESHDIR)/PatchData.cpp

MESHHEADERS = \
	$(patsubst $(MESHDIR)/%,$(localincludedir)/%,\
	$(wildcard $(MESHDIR)/*.hpp)) \
	$(patsubst $(MESHDIR)/%,$(localincludedir)/%,\
	$(wildcard $(MESHDIR)/*.hh))

MESHOBJ  = $(patsubst $(MESHDIR)/%.cpp, obj/%.o,$(MESHSRC))

ALL_MSQ_OBJ += $(MESHOBJ)
