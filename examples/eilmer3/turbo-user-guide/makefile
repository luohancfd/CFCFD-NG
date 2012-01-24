# makefile for turbo-user-guide
# PJ, Oct 2011

turbo-user-guide.pdf: turbo-user-guide.tex
	export BIBINPUTS=$HOME/papers/bibtex/
	pdflatex turbo-user-guide
	- bibtex turbo-user-guide
	- makeindex -g -l turbo-user-guide
	pdflatex turbo-user-guide
	pdflatex turbo-user-guide

clean:
	- rm turbo-user-guide.log \
		*.aux \
		*.backup \
		turbo-user-guide.idx \
		turbo-user-guide.ilg \
		turbo-user-guide.ind \
		turbo-user-guide.bbl \
		turbo-user-guide.blg \
		turbo-user-guide.out \
		turbo-user-guide.toc

all_clean:
	make clean
	- rm turbo-user-guide.pdf
