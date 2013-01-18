# makefile for estc report
# PJ, Jan 2013

estcj.pdf: estcj.tex
#	export BIBINPUTS=$HOME/papers/bibtex/
	pdflatex estcj
	- bibtex estcj
	pdflatex estcj
	pdflatex estcj

clean:
	- rm estcj.log estcj.aux estcj.bbl estcj.blg estcj.out \
		estcj.toc estcj.tex.backup estcj.tex~

all_clean:
	make clean
	- rm estcj.pdf
