
all:	pdf

illustrations.tex:	illustrations.Rnw
	~/bin/Sweave illustrations.Rnw

trees.tex:	trees.Rnw
	~/bin/Sweave trees.Rnw

clean:
	~/bin/texclean gpsi.tex
	rm -rf *pdf
	rm -rf *eps
	rm -rf illustrations.tex
	rm -rf trees.tex
	rm -rf multcomp.zip
	
pdf:	gpsi.tex illustrations.tex trees.tex references.bib header.tex
	texi2dvi --pdf gpsi.tex
	

zip:	
	~/bin/texclean gpsi.tex
	rm -rf multcomp.zip
	zip multcomp.zip *tex *pdf *bib