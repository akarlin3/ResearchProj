# Build the manuscript and supplement with tectonic.
# Figures are read from ../figures/manuscript/ (committed PDFs).
TECTONIC ?= tectonic

.PHONY: all manuscript supplement clean

all: manuscript supplement

manuscript: manuscript.tex
	$(TECTONIC) manuscript.tex

supplement: supplement.tex
	$(TECTONIC) supplement.tex

clean:
	rm -f manuscript.pdf supplement.pdf *.aux *.log *.out
