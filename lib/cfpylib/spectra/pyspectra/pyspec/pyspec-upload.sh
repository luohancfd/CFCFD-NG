#!/bin/bash
PYSPEC=pyspec
PYTHON=`which python`
DIST=pyspec-0.2
function upload_package {	

if [ -d "dist" ]; then
	echo "*** Removing dist directory ***"
	rm -rf dist
fi

$PYTHON setup.py build_sphinx upload_sphinx
$PYTHON setup.py sdist --formats=gztar,zip upload
scp dist/* stuwilkins,pyspec@frs.sourceforge.net:/home/frs/project/p/py/pyspec/$DIST

if [ -e "setup.cfg.inst" ]; then
	mv setup.cfg.inst setup.cfg
fi

}

function upload_docs {
	scp -r doc/redirect.html stuwilkins,pyspec@web.sourceforge.net:/home/groups/p/py/pyspec/htdocs/index.html
}

upload_package
upload_docs
