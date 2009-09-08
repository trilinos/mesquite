if ENABLE_DOXYGEN
  DOXY_SUBDIR = doxygen
endif
SUBDIRS = $(DOXY_SUBDIR)

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
	wrappers.tex \
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
	$(DVIPS) $< -o $@

users-guide.dvi: users-guide.tex $(UG_SRCS) 
	$(LATEX) $<
	$(BIBTEX) users-guide.aux
	$(LATEX) $<
	$(LATEX) $<

pdf: users-guide.pdf

ps: users-guide.ps


MOSTLYCLEANFILES = *.aux *.bbl *.blg *.lof *.log *.lot *.toc 
CLEANFILES = users-guide.dvi users-guide.pdf users-guide.ps


