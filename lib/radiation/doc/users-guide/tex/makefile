FILENAME = users-guide

$(FILENAME).pdf: $(FILENAME).tex \
		introduction.tex \
		installation.tex \
		getting-started.tex \
		input-files.tex \
		examples.tex
	pdflatex $(FILENAME).tex
	bibtex $(FILENAME)
	pdflatex $(FILENAME).tex
	pdflatex $(FILENAME).tex

clean:
	- rm *.aux *.bbl *.blg *.lof *.log *.lot *.toc *.out $(FILENAME).pdf
