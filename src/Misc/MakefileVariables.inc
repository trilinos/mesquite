MISCDIR = $(srcdir)/Misc

MISCSRC = $(MISCDIR)/MesquiteError.cpp \
	       $(MISCDIR)/MesquiteUtilities.cpp \
	       $(MISCDIR)/MsqTimer.cpp \
	       $(MISCDIR)/MsqMessage.cpp \
	       $(MISCDIR)/MsqHessian.cpp \
	       $(MISCDIR)/SimplifiedGeometryEngine.cpp \
##	       $(MISCDIR)/ExodusConversion.cpp

MISCHEADERS = \
	$(patsubst $(MISCDIR)/%,$(localincludedir)/%,\
	$(wildcard $(MISCDIR)/*.hpp)) \
	$(patsubst $(MISCDIR)/%,$(localincludedir)/%,\
	$(wildcard $(MISCDIR)/*.hh))

MISCOBJ  = $(patsubst $(MISCDIR)/%.cpp, obj/%.o, $(MISCSRC))

ALL_MSQ_OBJ += $(MISCOBJ)
