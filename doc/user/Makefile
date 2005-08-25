 # Commands
LATEX      = latex
DVI_TO_PDF = dvipdf
DVI_TO_PS  = dvips
BIBTEX     = bibtex

 # Input/output file basename
NAME = users-guide

 # aliases
pdf: $(NAME).pdf

ps:  $(NAME).ps

dvi: $(NAME).dvi

BIBLIOGRAPHY = mesquite.bib 

COMPONENTS = $(NAME).tex \
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
        trouble.tex 

$(NAME).pdf: $(NAME).dvi
	$(DVI_TO_PDF) $< $@

$(NAME).ps: $(NAME).dvi
	$(DVI_TO_PS) $< -o $@

$(NAME).dvi: $(COMPONENTS) $(BIBLIOGRAPHY)
	$(LATEX) $(NAME).tex
	$(BIBTEX) $(NAME).aux
	$(LATEX) $(NAME).tex
	$(LATEX) $(NAME).tex

TMPEXT = dvi lof log lot toc bbl
clean:
	rm -f $(COMPONENTS:%.tex=%.aux) $(TMPEXT:%=$(NAME).%)

dist-clean: clean
	rm -f $(NAME).ps $(NAME).pdf
