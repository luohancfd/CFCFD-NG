# makefile for generating source-code documentation

all: README-for-developers.html ./html/index.html
	echo "Should have finished building Developer notes."

clean:
	- rm -r ./html/* ./latex

./html/index.html:
	doxygen Doxyfile

README-for-developers.html: README-for-developers.txt
	asciidoc -b xhtml11 README-for-developers.txt
