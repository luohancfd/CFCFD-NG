#! /bin/bash
# prep-run-post.sh

e3prep.py --job=hydrogen --do-svg

e3shared.exe --job=hydrogen --run

e3post.py --job=hydrogen --output-file='eilmer3-computed.data' --slice-list=":,:,0,0"
