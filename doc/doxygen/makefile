# makefile for generating source-code documentation

all:
	doxygen Doxyfile

clean:
	- rm -r ./html/*.html ./html/*.png ./latex

tarfile :
	- cd $(HOME); mv cfd_doxy_doc.tar.gz cfd_doxy_doc.tar.gz.old
	cd $(HOME); tar --exclude ".svn" --exclude "*~" -cvf cfd_doxy_doc.tar cfcfd/code/doxygen
	cd $(HOME); gzip cfd_doxy_doc.tar
