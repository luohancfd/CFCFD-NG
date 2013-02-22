#!/bin/env python
# set_permissions.py
#
# Sets the permissions for the tunnel data files 
# for the entire tunnels collection.
#
# PJ, 28-May-2006
#
import os
import sys

def print_usage():
    print "Usage: set_permissions.py <top-directory>"
    return

if len(sys.argv) != 2:
    print_usage()
    sys.exit()

new_group_name = "tunnels"
top_dir = sys.argv[1]
print "set_permissions.py: Begin..."
print "   new_group_name:", new_group_name
print "   top_dir=", top_dir
file_count = 0
dir_count = 0
for (dir_path, dir_names, file_names) in os.walk(top_dir):
    print "dir_path=", dir_path
    os.chmod(dir_path, 0775)
    dir_count += 1
    for each_file in file_names:
        os.chmod(os.path.join(dir_path,each_file), 0664)
        file_count += 1
print "Done: ", file_count, "files and", dir_count, "directories."

