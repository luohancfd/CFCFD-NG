#!/bin/sh
# install.sh

INSTALL_DIR=${HOME}/e3bin
cd ..; rsync -av --exclude=".svn" --exclude="*~" cfpylib ${INSTALL_DIR}
