#!/bin/bash

gunzip test.dat.gz
onedval test.dat test.config one-d-props.txt
gzip test.dat

