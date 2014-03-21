#!/bin/bash

cd part1-inviscid
./run.sh
cd ..

cd part2-viscous
./run.sh
cd ..

cd part3-viscous-with-radiation
./run.sh
cd ..

echo "Done."
