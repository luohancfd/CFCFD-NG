#! /bin/sh
# file: setup_lut.sh

build-cea-lut --case=air

echo "We should now have a Look-Up-Table for air"

build-cea-lut --case=N2

echo "We should now have a Look-Up-Table for nitrogen"

build-cea-lut --case=air5species

echo "We should now have a Look-Up-Table for 5-species air"
