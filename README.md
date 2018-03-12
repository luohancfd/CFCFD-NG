# CFCFD-NG
This repository is forked from [CFCFD](http://cfcfd.mechmining.uq.edu.au/eilmer3.html). It only contains the version Eilmer3. The **master** branch is a direct clone of the [original repository](https://source.eqit.uq.edu.au/hg/cfcfd3) and **dev** branch contains my own modifications. I will try my best to keep both branches updated. Thanks again for the great work of original developers. 

## Changes:
### patch1
#### Tools
 - **ReadCellGrid.m:** Read Eilmer3 generated cell-centered grid. (You may not need it in your work)
 - **ReadGrid.m:** Read normal plot3d vertex-centered multiblock grid file (ext: g, xyz)
 - **ReadFlowEilmer3.m:** Read plot3d multiblock function file from Eilmer3 (ext: f, fmt)
 - **ReadFlow.m:** Read plot3d multiblock function file
 - **ReadName.m:** Read name file (ext: nam, name)
 - **ConcateBlocks.m:** Merge multiblock grid file and functionfile (Only test for my 2D cases which have uniformed divided multi-block grid)
 - **ExportXYZ.m/ExportFUNCTION.m:** Export correct plot3d format
 - **PLOT3D_post.m:** script to involke the post processing of previous sample
 - **poshax2tec.py**: Convert poshax3 output to tecplot format
#### Models
 - **VT relaxation model**: A new model called PolyFit, which can be used to for polynominal fit VT relaxation time or Eq. 5 in [paper](http://aip.scitation.org/doi/abs/10.1063/1.4813070). Check *energy-exchange-relaxation-time.cxx* for detail. A example input *o2-VT-QCT.lua* is included.
 - **Macheret-Fridman model**: Modify the models based on [paper](https://arc.aiaa.org/doi/abs/10.2514/1.T5375)
 
## Licenses:
> CFCFD program collection is a set of flow simulation tools for compressible fluids. Copyright (C) 1991-2017 Peter Jacobs, Rowan Gollan, Daniel Potter, Ingo Jahn, Anand Veeraragavan, Vince Wheatley, Daryl Bond, Chris James and other members of the CFCFD group.
> This collection is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.
> This program collection is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
> You should have received a copy of the [GNU-General-Public-License](https://www.gnu.org/licenses/gpl-3.0.html) along with this program. If not, see <http://www.gnu.org/licenses/>.
All the modifications I made is also licensed with [GPL 3.0](https://www.gnu.org/licenses/gpl-3.0.html). 
