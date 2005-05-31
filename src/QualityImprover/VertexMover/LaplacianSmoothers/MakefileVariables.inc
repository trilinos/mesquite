LAPLSMOOTHDIR = $(srcdir)/QualityImprover/VertexMover/LaplacianSmoothers

LAPLSMOOTHSRC = $(LAPLSMOOTHDIR)/LaplacianSmoother.cpp \
                $(LAPLSMOOTHDIR)/SmartLaplacianSmoother.cpp \
                $(LAPLSMOOTHDIR)/NormSquaredSmoother.cpp \
                $(LAPLSMOOTHDIR)/I_DFT_NoBarrierSmoother.cpp 

LAPLSMOOTHHEADERS = \
	$(patsubst $(LAPLSMOOTHDIR)/%,$(localincludedir)/%,\
	$(wildcard $(LAPLSMOOTHDIR)/*.hpp)) \
	$(patsubst $(LAPLSMOOTHDIR)/%,$(localincludedir)/%,\
	$(wildcard $(LAPLSMOOTHDIR)/*.hh))

LAPLSMOOTHOBJ  = $(patsubst $(LAPLSMOOTHDIR)/%.cpp, obj/%.o, $(LAPLSMOOTHSRC)) 

ALL_MSQ_OBJ += $(LAPLSMOOTHOBJ)
