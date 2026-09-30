# Build the design note. Independent of Docs/Makefile, which is the Sphinx tree.
MAIN = amrcomplex-design

.PHONY: all clean watch
all: $(MAIN).pdf

$(MAIN).pdf: $(wildcard *.tex) refs.bib
	latexmk -pdf -bibtex -halt-on-error -interaction=nonstopmode $(MAIN).tex

watch:
	latexmk -pdf -pvc $(MAIN).tex

clean:
	latexmk -C
	rm -f *.bbl *.blg
