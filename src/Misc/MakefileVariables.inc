MISCDIR = $(srcdir)/Misc

MISCSRC = $(MISCDIR)/MesquiteError.cpp \
	       $(MISCDIR)/MsqTimer.cpp \
	       $(MISCDIR)/MsqMessage.cpp \
	       $(MISCDIR)/MsqHessian.cpp \
	       $(MISCDIR)/PlanarDomain.cpp \
	       $(MISCDIR)/SphericalDomain.cpp

MISCHEADERS = \
	$(patsubst $(MISCDIR)/%,$(localincludedir)/%,\
	$(wildcard $(MISCDIR)/*.hpp)) \
	$(patsubst $(MISCDIR)/%,$(localincludedir)/%,\
	$(wildcard $(MISCDIR)/*.hh))

MISCOBJ  = $(patsubst $(MISCDIR)/%.cpp, obj/%.o, $(MISCSRC))

ALL_MSQ_OBJ += $(MISCOBJ)
