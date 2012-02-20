Flow (house-keeping) module
===========================

.. automodule:: cfpylib.flow

.. currentmodule:: cfpylib/flow/


blockflow2d
-----------

.. automodule:: cfpylib.flow.blockflow2d

.. autofunction:: apply_uflowz 

.. autoclass:: BlockFlow2D
   :members: __init__, init_arrays, read_solution_for_cell, read, data_in_dictionary


vtk_xml_writer
--------------

.. automodule:: cfpylib.flow.vtk_xml_writer

.. autofunction:: write_vtk_xml_unstructured_grid_2D

.. autofunction:: write_vtk_xml_multi_block_2D


tecplot_writer
--------------

.. automodule:: cfpylib.flow.tecplot_writer

.. autofunction:: write_tecplot_structured_block_2D

.. autofunction:: write_tecplot_file_per_block

.. autofunction:: write_tecplot_single_file


