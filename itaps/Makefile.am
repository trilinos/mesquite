BUILT_SOURCES += itaps/iMesh_protos.h \
                 itaps/iMeshP_protos.h \
                 itaps/iGeom_protos.h \
                 itaps/iRel_protos.h

EXTRA_DIST += $(top_srcdir)/itaps/mkprotos.sh \
              $(top_srcdir)/itaps/CMakeLists.txt

ITAPS_DIR = itaps
MSQ_INCLUDES += -I$(top_srcdir)/$(ITAPS_DIR)
INCLUDED_MAKEFILES += $(ITAPS_DIR)/Makefile.am

MSQ_HDRS += \
	$(ITAPS_DIR)/iBase.h \
	$(ITAPS_DIR)/iBase_FCDefs.h \
	$(ITAPS_DIR)/iMesh.h \
	$(ITAPS_DIR)/iMesh_protos.h \
	$(ITAPS_DIR)/iMeshP.h \
	$(ITAPS_DIR)/iMeshP_protos.h \
	$(ITAPS_DIR)/iGeom.h \
	$(ITAPS_DIR)/iGeom_protos.h \
	$(ITAPS_DIR)/iRel.h \
	$(ITAPS_DIR)/iRel_protos.h 

itaps/iMesh_protos.h: itaps/iMesh.h
	$(top_srcdir)/itaps/mkprotos.sh iMesh $< $@

itaps/iMeshP_protos.h: itaps/iMeshP.h
	$(top_srcdir)/itaps/mkprotos.sh iMeshP $< $@

itaps/iGeom_protos.h: itaps/iGeom.h
	$(top_srcdir)/itaps/mkprotos.sh iGeom $< $@

itaps/iRel_protos.h: itaps/iRel.h
	$(top_srcdir)/itaps/mkprotos.sh iRel $< $@


