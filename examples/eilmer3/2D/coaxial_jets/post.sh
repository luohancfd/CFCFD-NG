# Post-processing script for the coaxial_jets case by Wilson Y. K. Chan
# Modified by Sam Stennett, Nov 2014

job=coaxial_jets
tindx=last 

# vtk or visit files for visualisation
e3post.py --job=$job --tindx=$tindx --vtk-xml --add-mach --add-pitot-p

# Boundary layer slices to compare with experimental data
e3post.py --job=$job --tindx=$tindx --output-file=x000127.dat --slice-list="0,2,:,:;1,2,:,:;2,2,:,:;3,2,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x002.dat --slice-list="0,23,:,:;1,23,:,:;2,23,:,:;3,23,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x003.dat --slice-list="0,30,:,:;1,30,:,:;2,30,:,:;3,30,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x005.dat --slice-list="0,39,:,:;1,39,:,:;2,39,:,:;3,39,:,:"  --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x010.dat --slice-list="0,56,:,:;1,56,:,:;2,56,:,:;3,56,:,:"  --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x012.dat --slice-list="0,-2,:,:;1,-2,:,:;2,-2,:,:;3,-2,:,:"  --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x018.dat --slice-list="4,16,:,:;10,16,:,:;16,16,:,:;17,16,:,:;28,16,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x027.dat --slice-list="4,37,:,:;10,37,:,:;16,37,:,:;17,37,:,:;28,37,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x028.dat --slice-list="4,-1,:,:;10,-1,:,:;16,-1,:,:;17,-1,:,:;28,-1,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x042.dat --slice-list="5,25,:,:;11,25,:,:;18,25,:,:;19,25,:,:;29,25,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x043.dat --slice-list="5,26,:,:;11,26,:,:;18,26,:,:;19,26,:,:;29,26,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x062.dat --slice-list="6,13,:,:;12,13,:,:;20,13,:,:;21,13,:,:;30,13,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x081.dat --slice-list="6,34,:,:;12,34,:,:;20,34,:,:;21,34,:,:;30,34,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x082.dat --slice-list="6,35,:,:;12,35,:,:;20,35,:,:;21,35,:,:;30,35,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x101.dat --slice-list="7,13,:,:;13,13,:,:;22,13,:,:;23,13,:,:;31,13,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x102.dat --slice-list="7,14,:,:;13,14,:,:;22,14,:,:;23,14,:,:;31,14,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x121.dat --slice-list="7,29,:,:;13,29,:,:;22,29,:,:;23,29,:,:;31,29,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x123.dat --slice-list="7,31,:,:;13,31,:,:;22,31,:,:;23,31,:,:;31,31,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x151.dat --slice-list="8,11,:,:;14,11,:,:;24,11,:,:;25,11,:,:;32,11,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x153.dat --slice-list="8,12,:,:;14,12,:,:;24,12,:,:;25,12,:,:;32,12,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x181.dat --slice-list="8,30,:,:;14,30,:,:;24,30,:,:;25,30,:,:;32,30,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x190.dat --slice-list="8,37,:,:;14,37,:,:;24,37,:,:;25,37,:,:;32,37,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x220.dat --slice-list="9,14,:,:;15,14,:,:;26,14,:,:;27,14,:,:;33,14,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x258.dat --slice-list="9,36,:,:;15,36,:,:;26,36,:,:;27,36,:,:;33,36,:,:" --add-pitot-p --add-mach
e3post.py --job=$job --tindx=$tindx --output-file=x261.dat --slice-list="9,37,:,:;15,37,:,:;26,37,:,:;27,37,:,:;33,37,:,:" --add-pitot-p --add-mach
