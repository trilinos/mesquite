MESHDIR = $(srcdir)/Mesh

MESHSRCFILES = MeshSet.cpp \
               MsqMeshEntity.cpp \
               MsqVertex.cpp \
               PatchData.cpp \
               MeshImpl.cpp \
               TopologyInfo.cpp
               

MESHSRC = $(MESHSRCFILES:%=$(MESHDIR)/%)

MESHHEADERS = \
	$(patsubst $(MESHDIR)/%,$(localincludedir)/%,\
	$(wildcard $(MESHDIR)/*.hpp)) \
	$(patsubst $(MESHDIR)/%,$(localincludedir)/%,\
	$(wildcard $(MESHDIR)/*.hh))

MESHOBJ  = $(patsubst $(MESHDIR)/%.cpp, obj/%.o,$(MESHSRC))

ALL_MSQ_OBJ += $(MESHOBJ) $(TSTT_OBJS) 

TSTT_SRCS = $(MESHDIR)/MeshTSTT.cpp
