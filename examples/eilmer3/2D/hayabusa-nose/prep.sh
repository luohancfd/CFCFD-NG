#!/bin/bash
# prep.sh

cp ~/cfcfd2/lib/gas/reaction-schemes/air/TV-TE.lua .
cp ~/cfcfd2/lib/gas/reaction-schemes/air/Park93-s03-AIC-EIIC.lua .

e3prep.py --job=hayabusa --do-svg