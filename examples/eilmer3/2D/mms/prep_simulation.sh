#!/bin/bash
maxima --batch=make_source_terms.mac
python f90_to_lua.py
cp mms-regular.py mms.py
e3prep.py --job=mms

