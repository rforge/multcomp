
### Makefile for chapters of "A Handbook of Statistical Analyses using R"

BOOK = Multiple_Comparisons_in_R.tex
Rnw = $(shell ls *.Rnw)
tex = $(Rnw:.Rnw=.tex)

.SUFFIXES:
.SUFFIXES: .Rnw .tex

.Rnw.tex:
	echo "library(utils); Sweave(\"$<\")" | R --vanilla 

all: Echo texs book

book: $(tex) $(BOOK)
	texi2dvi --pdf $(BOOK)

Echo: ; @(echo Source files ; echo $(Rnw) ; echo $(tex) ; echo ''; echo $?. )
texs: $(tex)

clean:
	rm -f $(tex)
	rm -f *.pdf
	rm -f *~

texclean: clean
	rm -rf *.csv
	rm -rf *.rda
	rm -rf *.out
	rm -rf *.dvi
	rm -rf *.log
	rm -rf *.aux
	rm -rf *.bbl
	rm -rf *.blg
	rm -rf *.ind
	rm -rf *.idx
	rm -rf *.ilg
	rm -rf *.lof
	rm -rf *.lot
	rm -rf *.toc
	
archive: texclean
	rm multcomp_book.tgz
	tar --exclude=*svn* -czvf multcomp_book.tgz ../multcomp_book ../LaTeXBibTeX
	