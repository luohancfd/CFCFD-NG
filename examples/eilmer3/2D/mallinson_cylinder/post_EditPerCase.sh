# post.sh
# Usage: e3post.py [--help] [--job=<jobFileName>] [--tindx=<index|all>]
#                  [--vtk-xml] [--tecplot]
#                  [--ref-function=<python-script>]
#                  [--report-norms]
#                  [--compare-job=<jobFileName> [--compare-tindx=<index>]]
#                  [--output-file=<profile-data-file>]
#                  [--slice-list="blk-range,i-range,j-range,k-range;..."]
#                  [--slice-at-point="blk-range,index-pair,x,y,z;..."]
#                  [--surface-list="blk,surface-name;..."]
#                  [--add-pitot-p] [--add-total-p] [--add-mach]
#                  [--sample]


job=mallinson_cylinder
tindx=last #change as desired

# Produces vtk (Paraview & Visit readable) files
e3post.py --job=$job --tindx=$tindx --vtk-xml

# Produces slices of data from first adjacent cells along wall
# NEEDS TO BE MODIFIED DEPENDING ON NUMBER OF BLOCKS/CELLS IN CASE
#e3post.py --job=$job --tindx=$tindx --output-file="yline.dat" --slice-list="0,:,0,:;2,:,0,:;4,:,0,:;6,:,0,:;8,:,0,:;10,:,0,:;12,:,0,:;14,:,0,:"

# Produces boundary layer profile slice
# NEEDS TO BE MODIFIED DEPENDING ON NUMBER OF BLOCKS/CELLS IN CASE
#e3post.py --job=$job --tindx=$tindx --output-file="xline.dat" --slice-list="14,-18,:,:;15,-18,:,:"
