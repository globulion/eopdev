CREV = r2
PREV = r1
CNAME = paper-$(CREV)
PNAME = paper-$(PREV)
CSUPP = supplementary-$(CREV)
PSUPP = supplementary-$(PREV)

.DEFAULT: pdf bibtex

.PHONY: all clean pdf bibtex diff supp sdiff show ssupp

all: pdf bibtex show

exts := .aux .bbl .blg .lof .lot .log .pdf
clean:
	rm -rvf $(foreach ext,$(exts), $(CNAME)$(ext));
	rm -rvf $(foreach ext,$(exts), $(PNAME)$(ext));
	rm -rvf $(foreach ext,$(exts), $(CNAME)$(ext));
	rm -rvf $(foreach ext,$(exts), $(PNAME)$(ext));
	rm -rvf $(foreach ext,$(exts), $(CSUPP)$(ext));
	rm -rvf $(foreach ext,$(exts), $(PSUPP)$(ext));
	rm -rvf $(CNAME).diff.*
	rm -rvf $(CSUPP).diff.*

pdf:    $(CNAME).pdf
bibtex: $(CNAME).bbl
show: pdf
	acroread $(CNAME).pdf &

diff: pdf bibtex $(PNAME).pdf $(PNAME).bbl $(CNAME).diff.pdf

supp:  $(CSUPP).pdf $(CSUPP).bbl 
sdiff: supp $(PSUPP).pdf $(PSUPP).bbl $(CSUPP).diff.pdf

ssupp: supp
	acroread $(CSUPP).pdf &


$(CSUPP).diff.pdf: $(CSUPP).tex $(PSUPP).tex
	latexdiff --exclude-safecmd="bibitem,table" $(PNAME).tex $< > $(basename $<).diff.tex
	while (pdflatex $(basename $<).diff.tex; \
	  grep "Rerun to get cross" $(basename $<).diff.log ) do true ; \
	done
	while (bibtex $(basename $<).diff; pdflatex $(basename $<).diff.tex; \
	  grep "Rerun to get citations correct" $(basename $<).diff.log ) do true; \
	done
	bibtex $(basename $<).diff
	latexdiff --exclude-safecmd="bibitem" $(PNAME).bbl $(basename $<).bbl > $(basename $<).diff.bbl
	pdflatex $(basename $<).diff.tex
	bibtex $(basename $<).diff
	rm -rvf acs-$(basename $@).bib

$(CNAME).diff.pdf: $(CNAME).tex $(PNAME).tex
	latexdiff --exclude-safecmd="bibitem" $(PNAME).tex $< > $(basename $<).diff.tex
	while (pdflatex $(basename $<).diff.tex; \
	  grep "Rerun to get cross" $(basename $<).diff.log ) do true ; \
	done
	while (bibtex $(basename $<).diff; pdflatex $(basename $<).diff.tex; \
	  grep "Rerun to get citations correct" $(basename $<).diff.log ) do true; \
	done
	bibtex $(basename $<).diff
	latexdiff --exclude-safecmd="bibitem" $(PNAME).bbl $(basename $<).bbl > $(basename $<).diff.bbl
	pdflatex $(basename $<).diff.tex
	bibtex $(basename $<).diff
	rm -rvf acs-$(basename $@).bib

%.bbl: %.tex
	while (bibtex $(basename $<); pdflatex $<; \
	  grep "Rerun to get citations correct" $(basename $<).log ) do true; \
	done
	#rm -rvf acs-$(basename $@).bib

%.pdf: %.tex
	while (pdflatex $<; \
	  grep "Rerun to get cross" $(basename $<).log ) do true ; \
	done

