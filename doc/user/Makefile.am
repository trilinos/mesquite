SUBDIRS = doxygen

UG_SRCS = users-guide.tex \
	ack.tex \
	api.tex \
	caveats.tex \
	data.tex \
	extend.tex \
	intro.tex \
	mesh.tex \
	plan.tex \
	started.tex \
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
	$(DVIPS) $< -o $@

users-guide.dvi: $(UG_SRCS) 
	if "x$(srcdir)" = "x."; then \
	  $(LATEX) users-guide.tex ;\
	  $(BIBTEX) users-guide.aux ;\
	  $(LATEX) users-guide.tex ;\
	  $(LATEX) users-guide.tex ;\
	else \
	  bd=`pwd`; \
          (cd $(srcdir) ;\
	   $(LATEX) -output-directory=$$bd users-guide.tex ;\
	  $(BIBTEX) -output-directory=$$bd users-guide.aux ;\
	  $(LATEX) -output-directory=$$bd users-guide.tex ;\
	  $(LATEX) -output-directory=$$bd users-guide.tex \
	  ) \
	fi

pdf: users-guide.pdf

ps: users-guide.ps


MOSTLYCLEANFILES = *.aux *.bbl *.blg *.lof *.log *.lot *.toc
