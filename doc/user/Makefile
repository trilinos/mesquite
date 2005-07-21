 # Commands
LATEX      = latex
DVI_TO_PDF = dvipdf
DVI_TO_PS  = dvips

 # Input/output file basename
NAME = users-guide

 # aliases
pdf: $(NAME).pdf

ps:  $(NAME).ps

dvi: $(NAME).dvi


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
	$(DVI_TO_PDF) $< > $@

$(NAME).ps: $(NAME).dvi
	$(DVI_TO_PS) $< -o $@

$(NAME).dvi: $(NAME).toc $(COMPONENTS)
	$(LATEX) $(NAME).tex

# Need to invoke LaTeX twice for correct generation
# of table of contents in final document.  Invoke
# the first time here to generate the table of contents,
# removing the .dvi file generated without the correct
# table of contents.
$(NAME).toc: $(COMPONENTS)
	rm -f $@
	$(LATEX) $(NAME).tex
	rm -f $(NAME).dvi

TMPEXT = dvi lof log lot toc
clean:
	rm -f $(COMPONENTS:%.tex=%.aux) $(TMPEXT:%=$(NAME).%)
