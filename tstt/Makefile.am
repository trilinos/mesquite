TSTT_DIR = tstt
EXTRA_DIST += \
	$(TSTT_DIR)/TSTTB.sidl \
	$(TSTT_DIR)/TSTTG.sidl \
	$(TSTT_DIR)/TSTTM.sidl \
	$(TSTT_DIR)/TSTTR.sidl 

TSTT_BDIR = $(top_builddir)/tstt
if WITH_TSTT
  TSTT_INCLUDES = -I$(TSTT_BDIR)/base \
                  -I$(TSTT_BDIR)/mesh \
                  -I$(TSTT_BDIR)/geom \
                  -I$(TSTT_BDIR)/assoc 
  libmesquite_la_DEPENDENCIES += $(TSTT_HDRS)
  TSTT_HDRS = \
      $(TSTT_BDIR)/base/TSTTB.hh \
      $(TSTT_BDIR)/geom/TSTTG.hh \
      $(TSTT_BDIR)/mesh/TSTTM.hh \
      $(TSTT_BDIR)/assoc/TSTTR.hh 
else
  TSTT_INCLUDES =
  TSTT_HDRS =
endif
MSQ_INCLUDES += $(TSTT_INCLUDES)
BUILT_SOURCES += $(TSTT_HDRS)

$(TSTT_BDIR)/assoc/TSTTR.hh : $(TSTT_DIR)/TSTTR.sidl 
	$(BABEL) -cC++ -E -R$(TSTT_BDIR)/repo -o$(TSTT_BDIR)/assoc $<

$(TSTT_BDIR)/geom/TSTTG.hh : $(TSTT_DIR)/TSTTG.sidl 
	$(BABEL) -cC++ -E -R$(TSTT_BDIR)/repo -o$(TSTT_BDIR)/geom $<

$(TSTT_BDIR)/mesh/TSTTM.hh : $(TSTT_DIR)/TSTTM.sidl 
	$(BABEL) -cC++ -E -R$(TSTT_BDIR)/repo -o$(TSTT_BDIR)/mesh $<

$(TSTT_BDIR)/base/TSTTB.hh : $(TSTT_DIR)/TSTTB.sidl 
	$(BABEL) -tXML -o$(TSTT_BDIR)/repo $<
	$(BABEL) -cC++ -E -o$(TSTT_BDIR)/base $<

$(TSTT_BDIR)/assoc/TSTTR.hh : $(TSTT_BDIR)/geom/TSTTG.hh $(TSTT_BDIR)/mesh/TSTTM.hh
$(TSTT_BDIR)/geom/TSTTG.hh : $(TSTT_BDIR)/base/TSTTB.hh
$(TSTT_BDIR)/mesh/TSTTM.hh : $(TSTT_BDIR)/base/TSTTB.hh

clean-local:
	@for dir in repo assoc base geom mesh; do \
	  if test -d $(TSTT_BDIR)/$$dir; then rm -r $(TSTT_BDIR)/$$dir; fi; \
	done
