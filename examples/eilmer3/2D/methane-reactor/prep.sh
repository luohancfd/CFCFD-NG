#!/bin/bash
# prep.sh
echo "Methane combustion in a constant-volume reactor."

if [ -z ${CFCFD_SRC} ]; then CFCFD_SRC=${HOME}/cfcfd3; fi

INPUTFILES=${CFCFD_SRC}/lib/gas/reaction-schemes/methane-combustion
cp ${INPUTFILES}/grimech30.lua .
cp ${INPUTFILES}/grimech30.inp .
gasfile grimech30.inp thermally-perfect-grimech30.lua

echo "Run simple reactor calculation."
echo "This should take about half a minute."
./fixed_volume_reactor.py > cvfm-sp-output.dat

# cp ${CFCFD_SRC}/examples/perfectly-stirred-reactor/simple_reactor.py .
# python simple_reactor.py \
#     --type=cvfm \
#     --gas=thermally-perfect-grimech30.lua \
#     --chem=grimech30.lua \
#     --T=2000 \
#     --p=1 \
#     --X="{'CH4':1.0, 'O2':2.0, 'N2':7.52}" \
#     --tmax=4.0e-04 > cvfm-sp-output.dat

e3prep.py --job=psr
exit $?
