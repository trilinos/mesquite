 # Commands
LATEX      = latex
DVI_TO_PDF = dvipdf
DVI_TO_PS  = dvips

 # Input/output file basename
NAME = users-guide

 # aliases
pdf: $(NAME).pdf

ps:  $(MAME).ps

dvi: $(NAME).dvi



$(NAME).pdf: $(NAME).dvi
	$(DVI_TO_PDF) $< $@

$(NAME).ps:
	$(DVI_TO_PS) $< $@

$(NAME).dvi: $(NAME).toc $(NAME).tex
	$(LATEX) $(NAME).tex

# Need to invoke LaTeX twice for correct generation
# of table of contents in final document.  Invoke
# the first time here to generate the table of contents,
# removing the .dvi file generated without the correct
# table of contents.
$(NAME).toc: $(NAME).tex
	rm -f $@
	$(LATEX) $<
	rm -f $(NAME).dvi

TMPEXT = dvi lof log lot toc
clean:
	rm -f *.aux $(TMPEXT:%=$(NAME).%)
