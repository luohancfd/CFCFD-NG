#!/usr/bin/env python
#
# author: Rowan J. Gollan
# date: 15-Jun-2007

import sys

try:
    from xml.etree.ElementTree import ElementTree, Element, SubElement
except ImportError:
    from elementtree.ElementTree import ElementTree, Element, SubElement
except:
    print "It is not possible to use module: vtk_xml_writer."
    print "without ElementTree package installed as part of your Python interpreter."
    print "Module not loaded."
    raise ImportError("No ElementTree module available.")


def write_vtk_xml_unstructured_grid_2D(grid, block, filename, data_list):
    """
    This function writes a serial vtk file in XML format for an unstructured grid.

    Inputs
    ------
    grid      : A BlockGrid2D object
    block     : A BlockFlow2D onject which corresponds to grid
    filename  : A filename for the output file (without extension)
    data_list : A list of keys for the desired data from the block

    """
    num_points = grid.ni * grid.nj
    num_cells = block.ni * block.nj
        
    vtkfile = Element("VTKFile", type="UnstructuredGrid")
    u_grid = SubElement(vtkfile, "UnstructuredGrid")
    piece = SubElement(u_grid, "Piece", NumberOfPoints=str(num_points), NumberOfCells=str(num_cells))
    points = SubElement(piece, "Points")
    d_array = SubElement(points, "DataArray", type="Float32", NumberOfComponents="3", format="ascii")
    d_array.text = ""
    for i in range(grid.ni):
        for j in range(grid.nj):
            d_array.text += "%20.12e %20.12e %20.12e\n" % (grid.x[i][j], grid.y[i][j], 0.0)

    cells = SubElement(piece, "Cells")
    d_array = SubElement(cells, "DataArray", type="Int32", Name="connectivity")
    d_array.text = ""
    for i in range(block.ni):
        for j in range(block.nj):
            d_array.text += "%d %d %d %d \n" % (i*grid.nj+j, (i+1)*grid.nj+j,
                                                (i+1)*grid.nj+(j+1), i*grid.nj+(j+1))

    d_array = SubElement(cells, "DataArray", type="Int32", Name="offsets")
    d_array.text = ""
    for i in range(block.ni):
        for j in range(block.nj):
            d_array.text += "%d \n" % ( 4*(i*block.nj + j + 1) )

    d_array = SubElement(cells, "DataArray", type="UInt8", Name="types")
    d_array.text = ""
    for i in range(block.ni):
        for j in range(block.nj):
            d_array.text += "%d \n" % (9) # VTK_quad == 9

    data = block.data_in_dictionary()
    s_labels = ""
    v_labels = ""
    for item, label, d_type in data_list:
        if d_type == "s":
            s_labels += label + " "
        elif d_type == "v":
            v_labels += label + " "
        else:
            print "Unknown data type:", d_type
            print "Expected:"
            print " s - scalar quantity"
            print " v - vector quantity"

    if s_labels != "" and v_labels != "":
        # Both scalar and vector quantities present
        cell_data = SubElement(piece, "CellData", Scalars=s_labels, Vectors=v_labels)
    elif s_labels != "" and v_labels == "":
        # Only scalar data present
        cell_data = SubElement(piece, "CellData", Scalars=s_labels)
    elif v_labels != "" and s_labels == "":
        # Only vector data present
        cell_data = SubElement(piece, "CellData", Vectors=v_labels)
    else:
        # No data requested
        print "No data written as no data was requested."
        return
    
    for item, label, d_type in data_list:
        if d_type == "s":
            d_array = SubElement(cell_data, "DataArray", type="Float32", Name=label)
            d_array.text = ""
            for i in range(block.ni):
                for j in range(block.nj):
                    d_array.text += "%20.12e \n" % (data[item][i][j])
        elif d_type == "v":
            d_array = SubElement(cell_data, "DataArray", type="Float32", Name=label,
                                 NumberOfComponents="3")
            d_array.text = ""
            for i in range(block.ni):
                for j in range(block.nj):
                    d_array.text += "%20.12e %20.12e %20.12e \n" % (
                        data[item][0][i][j], data[item][1][i][j], data[item][2][i][j])
        else:
            pass

    ElementTree(vtkfile).write("%s.%s" % (filename, "vtu"))
    return

def write_vtk_xml_multi_block_2D( grid_and_blocks, filename, data_list):
    """
    This function writes a multi block of flow data as individual files
    and a coordinating file.

    Inputs
    ------
    grid_and_blocks  : A list of tuples containing (grid, block) objects.
    filename         : The base name for the collection of files.
    data_list        : A list of tuples of the form (symbol, name)
    """

    for i, (grid, block) in enumerate(grid_and_blocks):
        fname = filename + "_%03d" % (i)
        write_vtk_xml_unstructured_grid_2D(grid, block, fname, data_list)

    vtkfile = Element("VTKFile", type="PUnstructuredGrid")
    pu_grid = SubElement(vtkfile, "PUnstructuredGrid", GhostLevel="0")
    ppoints = SubElement(pu_grid, "PPoints")
    pd_array = SubElement(ppoints, "PDataArray", type="Float32", NumberOfComponents="3")
    
    s_labels = ""
    v_labels = ""
    for item, label, d_type in data_list:
        if d_type == "s":
            s_labels += label + " "
        elif d_type == "v":
            v_labels += label + " "
        else:
            print "Unknown data type:", d_type
            print "Expected:"
            print " s - scalar quantity"
            print " v - vector quantity"

    if s_labels != "" and v_labels != "":
        # Both scalar and vector quantities present
        pcell_data = SubElement(pu_grid, "PCellData", Scalars=s_labels, Vectors=v_labels)
    elif s_labels != "" and v_labels == "":
        # Only scalar data present
        pcell_data = SubElement(pu_grid, "PCellData", Scalars=s_labels)
    elif v_labels != "" and s_labels == "":
        # Only vector data present
        pcell_data = SubElement(pu_grid, "PCellData", Vectors=v_labels)
    else:
        # No data requested
        print "No data written as no data was requested."
        return

    for item, label, d_type in data_list:
        if d_type == "s":
            pd_array = SubElement(pcell_data, "PDataArray", type="Float32", Name=label)
        elif d_type == "v":
            pd_array = SubElement(pcell_data, "PDataArray", type="Float32", Name=label,
                                  NumberOfComponents="3")
        else:
            pass
        
    for i in range(len(grid_and_blocks)):
        fname = filename + "_%03d.vtu" % (i)
        piece = SubElement(pu_grid, "Piece", Source=fname)

    ElementTree(vtkfile).write("%s.%s" % (filename, "pvtu"))
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
    
    write_vtk_xml_multi_block_2D( mblocks, "cone20_test", [("p", "pressure", "s"),
                                                           ("T", "temperature", "s"),
                                                           ("vel", "velocity", "v")])
    print "Done."

    
    
    
    
        
    
    
