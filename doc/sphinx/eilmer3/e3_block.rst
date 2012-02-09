e3_block
========

.. automodule:: e3_block

Generic functions and classes
-----------------------------

.. autofunction:: identify_block_connections

.. autofunction:: e3_block.make_patch

.. autofunction:: e3_block.close_enough

.. autoclass:: e3_block.Block

.. automethod:: e3_block.Block.set_BC

.. automethod:: e3_block.Block.set_WBC

Two-dimensional flows
---------------------

.. autoclass:: e3_block.Block2D

.. automethod:: e3_block.Block2D.__init__

.. automethod:: e3_block.Block2D.cell_centre_location

.. autofunction:: connect_blocks_2D

.. autofunction:: identify_block_connections_2D

.. autoclass:: e3_block.MultiBlock2D

.. automethod:: e3_block.MultiBlock2D.__init__

.. autoclass:: e3_block.SuperBlock2D

.. automethod:: e3_block.SuperBlock2D.__init__

Three-dimensional flows
-----------------------

.. autoclass:: e3_block.Block3D

.. automethod:: e3_block.Block3D.__init__

.. automethod:: e3_block.Block3D.cell_centre_location

.. autoclass:: e3_block.MultiBlock3D

.. automethod:: e3_block.MultiBlock3D.__init__

.. autoclass:: e3_block.SuperBlock3D

.. automethod:: e3_block.SuperBlock3D.__init__

.. autofunction:: identify_block_connections_3D

.. autofunction:: identify_colocated_vertices

.. autofunction:: connect_blocks_3D

.. autofunction:: cell_count_consistent_3D

