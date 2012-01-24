#!/bin/sh
# install.sh

INSTALL_DIR=${HOME}/cfd_bin
cd ..; rsync -av --exclude=".svn" --exclude="*~" cfpylib ${INSTALL_DIR}
