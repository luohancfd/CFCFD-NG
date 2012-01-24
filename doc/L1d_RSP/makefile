# makefile for Riemann sub-problem report
# DFP, Dec 2006

rsp.pdf : rsp.dvi
	dvips -o rsp.ps rsp
	ps2pdf -sPAPERSIZE=a4 rsp.ps rsp.pdf
	- rm rsp.ps.gz
	gzip rsp.ps

rsp.dvi: rsp.tex
	latex rsp
	- bibtex rsp
	latex rsp
	latex rsp

clean:
	- rm rsp.log \
		*.aux \
		rsp.bbl \
		rsp.blg \
		rsp.toc \
		rsp.dvi

all_clean:
	make clean
	- rm rsp.ps.gz rsp.pdf
