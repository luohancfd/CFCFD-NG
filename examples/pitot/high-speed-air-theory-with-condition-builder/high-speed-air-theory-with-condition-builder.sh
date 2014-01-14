#! /bin/sh
# file: high-speed-air-theory-with-condition-builder.sh
# Shell script to run the fully theoretical 
# air example for pitot that also uses the condition building function.
# This is a short example (24 tests), but a large one will take a while to run
#Chris James (c.james4@uq.edu.au) - 14/01/14

pitot_condition_builder.py --config_file=high-speed-air-theory-with-condition-builder.cfg

echo "Done"

