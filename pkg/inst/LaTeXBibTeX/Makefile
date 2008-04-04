
### Makefile for BibTeX file of "Multiple Comparisons in R"

MainSrc = $(shell ls  *.in)
Bib = $(MainSrc:.in=.bib)

.SUFFIXES:
.SUFFIXES: .in .bib

.in.bib:
	echo "library(utils); library("HSAUR"); HSAUR:::pkgversions(\"$<\")" | R --vanilla > $@

all: Echo Bibs

Echo: ; @(echo Source files ; echo $(MainSrc) ; echo '')
Bibs: $(Bib)

clean:
	rm -f $(Bib)
	rm -f *~
