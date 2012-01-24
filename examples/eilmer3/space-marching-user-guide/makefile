# makefile for eilmer3-k-omega-user-guide
# WYKC & PJ, 2011

eilmer3-user-guide.dvi: eilmer3-space-marching.tex \
                       ./chap1-introduction/introduction.tex \
                       ./chap2-stream-tube/stream-tube.tex \
                       ./chap3-axi-scramjet/axi-scramjet.tex
#	export BIBINPUTS=$HOME/papers/bibtex/
	pdflatex eilmer3-space-marching
	- bibtex eilmer3-space-marching
	pdflatex eilmer3-space-marching
	pdflatex eilmer3-space-marching

clean:
	- rm eilmer3-space-marching.log \
		*.aux \
		eilmer3-space-marching.bbl \
		eilmer3-space-marching.blg \
		eilmer3-space-marching.toc \
		./chap1-introduction/introduction.aux \
		./chap2-stream-tube/stream-tube.aux \
		./chap3-axi-scramjet/axi-scramjet.aux

all_clean:
	make clean
	- rm eilmer3-space-marching.pdf
