# ss3_post.sh
# By default, e3post.py grabs the solution at final time.
e3post.py --job=ss3 --vtk-xml

e3post.py --job=ss3 --slice-list="0,:,0,:" --output-file=ss3_stag_line.data
awk -f locate_shock.awk ss3_stag_line.data > ss3.result

