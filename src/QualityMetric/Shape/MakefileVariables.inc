QMSHAPEDIR = $(srcdir)/QualityMetric/Shape

QMSHAPESRC = $(QMSHAPEDIR)/AspectRatioGammaQualityMetric.cpp \
             $(QMSHAPEDIR)/ConditionNumberQualityMetric.cpp \
             $(QMSHAPEDIR)/GeneralizedConditionNumberQualityMetric.cpp \
	     $(QMSHAPEDIR)/MeanRatioQualityMetric.cpp 

QMSHAPEHEADERS = \
	$(patsubst $(QMSHAPEDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QMSHAPEDIR)/*.hpp)) \
	$(patsubst $(QMSHAPEDIR)/%,$(localincludedir)/%,\
	$(wildcard $(QMSHAPEDIR)/*.hh))

QMSHAPEOBJ  = $(patsubst $(QMSHAPEDIR)/%.cpp, obj/%.o, $(QMSHAPESRC)) 
