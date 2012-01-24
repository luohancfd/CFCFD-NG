# makefile for eilmer3-k-omega-user-guide
# WYKC & PJ, 2011

eilmer3-user-guide.dvi: eilmer3-turbulence-model-validation.tex \
                       ./chap1-introduction/introduction.tex \
                       ./chap2-flat-plate/flat-plate.tex \
                       ./chap3-mallinson-cylinder/mallinson-cylinder.tex \
                       ./chap4-backward-facing-step/backward-facing-step.tex \
                       ./chap5-coaxial-jets/coaxial-jets.tex \
                       ./chap6-conclusion/conclusion.tex \
                       ./evaluation-template/evaluation-template.tex 
#	export BIBINPUTS=$HOME/papers/bibtex/
	pdflatex eilmer3-turbulence-model-validation
	- bibtex eilmer3-turbulence-model-validation
	pdflatex eilmer3-turbulence-model-validation
	pdflatex eilmer3-turbulence-model-validation

clean:
	- rm eilmer3-turbulence-model-validation.log \
		*.aux \
		eilmer3-turbulence-model-validation.bbl \
		eilmer3-turbulence-model-validation.blg \
		eilmer3-turbulence-model-validation.toc \
                ./chap1-introduction/introduction.aux \
                ./chap2-flat-plate/flat-plate.aux \
                ./chap3-mallinson-cylinder/mallinson-cylinder.aux \
                ./chap4-backward-facing-step/backward-facing-step.aux \
                ./chap5-coaxial-jets/coaxial-jets.aux \
                ./chap6-conclusion/conclusion.aux \
                ./evaluation-template/evaluation-template.aux

all_clean:
	make clean
	- rm eilmer3-turbulence-model-validation.pdf
