#!/usr/bin/env python
"""
tecplot_writer.py: Writing of BlockGrid2D and BlockFlow2D (mbcns2) objects to VTK XML files.

Author: Rowan Gollan

Date: 22-Oct-2007
"""

import sys

def write_tecplot_structured_block_2D(grid, block, filename, data_list):
    """
    This function writes an unstructured grid in tecplot format.
    
    :param grid: A BlockGrid2D object
    :param block: A BlockFlow2D object which corresponds to grid
    :param filename: A filename for the output file (without extension)
    :data_list: A list of keys for the desired data from the block
    """

    # Write header information
    f = open(filename + ".tec", 'a')
    f.write("TITLE = \"%s\"\n" % grid.label)
    f.write("VARIABLES = \"X\", \"Y\"")
    for i, (item, label, d_type) in enumerate(data_list):
        if d_type == "s":
            f.write(", \"%s\"" % (label))
        elif d_type == "v":
            print "Tecplot does not have a vector data_type."
            print "Plotting %s data using scalars instead." % (label)

            data_list.append((item+"_x", item+"_x", "s"))
            data_list.append((item+"_y", item+"_y", "s"))
            data_list.append((item+"_z", item+"_z", "s"))
            #f.write(", \"%s\", \"%s\", \"%s\"" % (item+"_x", item+"_y", item+"_z"))
        else:
            print "Unknown data type:", d_type
            print "Expected:"
            print " s - scalar quantity"
            print " v - vector quantity"
    f.write("\n")
    
    # remove vector types
    for i, (item, label, d_type) in enumerate(data_list):
        if d_type == "v": data_list.remove((item, label, d_type))

    f.write("ZONE I=%i, J=%i, DATAPACKING=BLOCK\n" % (grid.ni, grid.nj))

    data_loc = ""
    for i, (item, label, d_type) in enumerate(data_list):
        data_loc += (str(3+i)+"=CELLCENTERED, ")
    if data_loc != "":
        f.write("VARLOCATION=(%s)\n" % (data_loc))

    # Write grid information
    for j in range(grid.nj):
        for i in range(grid.ni):
            f.write("%e " % (grid.x[i][j]))
        f.write("\n")
        
    for j in range(grid.nj):
        for i in range(grid.ni):
            f.write("%e " %(grid.y[i][j]))
        f.write("\n")

    # Write data information
    data = block.data_in_dictionary()
    for item, label, d_type in data_list:
        for j in range(block.nj):
            for i in range(block.ni):
                f.write("%20.12e " % (data[item][i][j]))
            f.write("\n")
    
    f.close()

    #print "grid in %s should be %d lines" % (filename, grid.nj*2.0)
    #print "data in %s should be %d lines of %d values long" % (filename, len(data_list), block.nj*block.ni)

    return

def write_tecplot_file_per_block(grid_and_blocks, filename, data_list):
    """
    This function writes a structured grid in tecplot format.
    
    :param grid_and_blocks: A list of tuples containing (grid, block) objects.
    :param filename: A filename for the output file (without extension)
    :param data_list: A list of keys for the desired data from the block
    """

    for i, (grid, block) in enumerate(grid_and_blocks):
        fname = filename + "_%03d" % (i)
        write_tecplot_structured_block_2D(grid, block, fname, data_list)
    
    return

def write_tecplot_single_file(grid_and_blocks, filename, data_list):
    """
    This function writes a structured grid in tecplot format.
    
    :param grid_and_blocks: A list of tuples containing (grid, block) objects.
    :param filename: A filename for the output file (without extension)
    :param data_list: A list of keys for the desired data from the block
    """

    for i, (grid, block) in enumerate(grid_and_blocks):
        fname = filename
        write_tecplot_structured_block_2D(grid, block, fname, data_list)
    
    return

if __name__ == '__main__':
    from cfpylib.flow.blockflow2d import BlockFlow2D
    from blockgrid2d import BlockGrid2D
    import sys
    print "Begin read block data for cone20 test case..."
    try:
        file_pointer = open("cone20.s", "r")
    except:
        print "This test only works when copied to the cone20 directory."
        print "~/cfcfd2/examples/mbcns2/cone20"
        sys.exit(1)
    block0 = BlockFlow2D()
    block1 = BlockFlow2D()
    for it in range(4):
        block0.read(file_pointer)
        block1.read(file_pointer)
        print "Have read data for time t=", block0.t

    print "Begin read grids..."
    gfp = open("cone20.g", "r")

    grid0 = BlockGrid2D()
    grid1 = BlockGrid2D()
    grid0.read(gfp)
    grid1.read(gfp)

    mblocks = [(grid0, block0), (grid1, block1)]
    data_list = [("p", "pressure", "s"), ("T", "temperature", "s"), ("vel", "velocity", "v")]
    
    write_tecplot_single_file(mblocks, "cone20_test", data_list)

    print "Done."
                                 
