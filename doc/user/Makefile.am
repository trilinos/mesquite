SUBDIRS = doxygen

UG_SRCS = users-guide.tex \
	ack.tex \
	api.tex \
	caveats.tex \
	data.tex \
	geom.tex \
	extend.tex \
	intro.tex \
	mesh.tex \
	plan.tex \
	install.tex \
	examples.tex \
	support.tex \
	team.tex \
        trouble.tex \
        MesquiteUI.eps \
        hole_in_square.ps \
        msq-paradigm.eps \
        square_rand.eps \
        mesquite.bib

EXTRA_DIST = $(UG_SRCS)

if WITH_USER_GUIDE_PDF
  pdf_DATA = users-guide.pdf
endif
if WITH_USER_GUIDE_PS
  ps_DATA = users-guide.ps
endif

SUFFIXES = .dvi .pdf .ps
.ps.pdf:
	$(PS2PDF) $< $@

.dvi.ps:
	bd=`pwd`; (cd $(srcdir) && $(DVIPS) $$bd/$< -o $$bd/$@)

users-guide.dvi: $(UG_SRCS) 
	for file in $(UG_SRCS); do if test ! -f $$file -o $(srcdir)/$$file -nt ./$$file; then cp $(srcdir)/$$file .; fi; done
	($(LATEX) users-guide.tex && \
	 $(BIBTEX) users-guide.aux && \
	 $(LATEX) users-guide.tex && \
	 $(LATEX) users-guide.tex); \
        exitval=$$? ;\
        if test ! $(srcdir) -ef .; then rm $(UG_SRCS); fi ;\
	exit $$exitval

pdf: users-guide.pdf

ps: users-guide.ps


MOSTLYCLEANFILES = *.aux *.bbl *.blg *.lof *.log *.lot *.toc 
CLEANFILES = users-guide.dvi users-guide.pdf users-guide.ps


