FILENAME = users-guide

$(FILENAME).pdf: $(FILENAME).tex \
		introduction.tex \
		thermodynamics.tex \
		trans-prop.tex \
		chem-kinetics.tex \
		installation.tex \
		getting-started.tex \
		input-files.tex \
		cfd.tex \
		python.tex \
		examples.tex \
		api.tex \
		lut.tex \
		user-defined.tex \
		poly-curves.tex \
		code.tex
	pdflatex $(FILENAME).tex
	bibtex $(FILENAME)
	pdflatex $(FILENAME).tex
	pdflatex $(FILENAME).tex

clean:
	- rm *.aux *.bbl *.blg *.lof *.log *.lot *.toc *.out $(FILENAME).pdf
