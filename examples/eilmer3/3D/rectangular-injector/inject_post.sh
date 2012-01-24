#! /bin/sh
# inject_post.sh
e3post.py --job=inject --vtk-xml
python mix_and_vort.py inject mixvort.dat
