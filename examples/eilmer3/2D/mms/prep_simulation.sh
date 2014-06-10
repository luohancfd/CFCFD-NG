#!/bin/bash
python make_source_terms.py
cp mms-regular.py mms.py
e3prep.py --job=mms

