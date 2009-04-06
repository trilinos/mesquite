BUILT_SOURCES += itaps/iMesh_protos.h \
                 itaps/iGeom_protos.h \
                 itaps/iRel_protos.h

EXTRA_DIST = mkprotos.sh

itaps/iMesh_protos.h: itaps/iMesh.h
	$(top_srcdir)/itaps/mkprotos.sh iMesh $< $@

itaps/iGeom_protos.h: itaps/iGeom.h
	$(top_srcdir)/itaps/mkprotos.sh iGeom $< $@

itaps/iRel_protos.h: itaps/iRel.h
	$(top_srcdir)/itaps/mkprotos.sh iRel $< $@
