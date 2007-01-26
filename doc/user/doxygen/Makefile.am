SRC_FILE = Mesquite-user.dox
DOXYGEN_HTML = @DOXYGEN_HTML@
DOXYGEN_PDF = @DOXYGEN_PDF@
DOXYGEN_PS = @DOXYGEN_PS@
USER_GUIDE_PDF = @USER_GUIDE_PDF@
USER_GUIDE_PS = @USER_GUIDE_PS@

if WITH_DOXYGEN_PDF
  DOXYGEN_LATEX=yes
  pdf_DATA = MesquiteAPI.pdf
else
  pdf_DATA =
if WITH_DOXYGEN_PS
    DOXYGEN_LATEX=yes
else
    DOXYGEN_LATEX=no
endif
endif

if WITH_DOXYGEN_PS
  ps_DATA = MesquiteAPI.ps
else
  ps_DATA =
endif

if WITH_DOXYGEN_HTML
  html_DATA = html/*
else
  html_DATA =
endif

guide.dox: $(SRC_FILE) Makefile
	$(SED)	-e 's/^GENERATE_HTML *=.*$$/GENERATE_HTML = $(DOXYGEN_HTML)/' \
		-e 's/^GENERATE_HTMLHELP *=.*$$/GENERATE_HTMLHELP = no/' \
		-e 's/^GENERATE_CHI *=.*$$/GENERATE_CHI = no/' \
		-e 's/^GENERATE_LATEX *=.*$$/GENERATE_LATEX = $(DOXYGEN_LATEX)/' \
		-e 's/^GENERATE_RTF *=.*$$/GENERATE_RTF = no/' \
		-e 's/^GENERATE_MAN *=.*$$/GENERATE_MAN = no/' \
		-e 's/^GENERATE_XML *=.*$$/GENERATE_XML = no/' \
		$< > $@
	echo "INPUT = $(top_srcdir)/src $(top_srcdir)/include $(top_builddir)/include" >> $@



html/index.html latex/refman.tex: guide.dox
	$(DOXYGEN) $<

latex/refman.dvi: latex/refman.tex latex/doxygen.sty
	cd latex && $(LATEX) refman.tex	
	cd latex && makeindex refman.tex
	cd latex && $(LATEX) refman.tex	
	cd latex && $(LATEX) refman.tex	
	cd latex && $(LATEX) refman.tex	
	cd latex && $(LATEX) refman.tex	

latex/refman.ps: latex/refman.dvi
	cd latex && $(DVIPS) refman.dvi

MesquiteAPI.ps: latex/refman.ps
	cp $< $@

MesquiteAPI.pdf: MesquiteAPI.ps
	ps2pdf $< $@

