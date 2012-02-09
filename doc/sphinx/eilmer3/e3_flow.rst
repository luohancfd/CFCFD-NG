e3_flow
=======

.. automodule:: e3_flow

Preparing flow field files
--------------------------

FlowCondition
^^^^^^^^^^^^^

.. autoclass:: e3_flow.FlowCondition

.. automethod:: e3_flow.FlowCondition.__init__

.. automethod:: e3_flow.FlowCondition.to_dict

.. automethod:: e3_flow.FlowCondition.write_to_ini_file

.. autofunction:: variable_list_for_cell

.. autofunction:: write_cell_data

StructuredGridFlow
^^^^^^^^^^^^^^^^^^

.. autoclass:: e3_flow.StructuredGridFlow

.. automethod:: e3_flow.StructuredGridFlow.read

.. automethod:: e3_flow.StructuredGridFlow.get_cell_data

.. automethod:: e3_flow.StructuredGridFlow.find_nearest_cell_centre

.. automethod:: e3_flow.StructuredGridFlow.add_aux_variables

.. automethod:: e3_flow.StructuredGridFlow.write_gnuplot_header

.. automethod:: e3_flow.StructuredGridFlow.write_gnuplot_data_for_cell

House-keeping
^^^^^^^^^^^^^

.. autofunction:: read_all_blocks

.. autofunction:: add_auxiliary_variables

.. autofunction:: locate_cell_and_block

ExistingSolution
^^^^^^^^^^^^^^^^

.. autoclass:: e3_flow.ExistingSolution

.. automethod:: e3_flow.ExistingSolution.__init__

.. automethod:: e3_flow.ExistingSolution.interpolate_flow_condition


Writing plot files
------------------

.. autofunction:: uflowz

.. autofunction:: write_VTK_XML_unstructured_file

.. autofunction:: write_VTK_XML_files

.. autofunction:: write_plot3d_files

Profile extraction and writing
------------------------------

.. autofunction:: decode_range_from_string

.. autofunction:: write_profile_data

.. autofunction:: convert_string

.. autofunction:: write_profile_along_line

Comparing with a reference function
-----------------------------------

.. autofunction:: compute_difference_in_flow_data

.. autofunction:: compute_difference_in_flow_data2

.. autofunction:: compute_volume_weighted_norms

.. autofunction:: pretty_print_norms

Radiation post-processing
-------------------------

.. autofunction:: tangent_slab_along_slice
