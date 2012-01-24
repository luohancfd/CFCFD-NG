# ss3_run_py.sh
# Shell script to set up and run Sawada & Dendou's sphere case 3.

# For a clean start
e3prep.py --job=ss3.py --do-svg

# The main event
time e3shared.exe --job=ss3 --run


