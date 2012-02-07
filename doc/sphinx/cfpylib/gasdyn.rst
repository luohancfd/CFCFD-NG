Gas Dynamics module
===================

.. automodule:: cfpylib.gasdyn

.. currentmodule:: cfpylib/gasdyn

billig
------

.. automodule:: cfpylib.gasdyn.billig

.. autofunction:: delta_over_R

.. autofunction:: Rc_over_R

.. autofunction:: x_from_y

.. autofunction:: y_from_x


cea2_gas
--------

.. automodule:: cfpylib.gasdyn.cea2_gas

.. autofunction:: test_for_cea_exe

.. autofunction:: run_cea_program

.. autofunction:: get_cea2_float

.. autoclass:: Gas
   :members: __init__, get_eq_massf, get_eq_molef, set_from_pAndT, write_state, EOS, Shock

ideal_gas_flow
--------------

.. automodule:: cfpylib.gasdyn.ideal_gas_flow

.. autofunction:: A_Astar

.. autofunction:: T0_T

.. autofunction:: p0_p

.. autofunction:: r0_r

.. autofunction:: m2_shock

.. autofunction:: r2_r1

.. autofunction:: u2_u1

.. autofunction:: p2_p1

.. autofunction::T2_T1 

.. autofunction:: p02_p01

.. autofunction:: DS_Cv

.. autofunction:: T0_T0star

.. autofunction:: M_Rayleigh

.. autofunction:: T_Tstar

.. autofunction:: p_pstar

.. autofunction:: r_rstar

.. autofunction:: p0_p0star

.. autofunction:: PM1

.. autofunction:: PM2

.. autofunction:: beta_obl

.. autofunction:: theta_obl

.. autofunction:: M2_obl

.. autofunction:: r2_r1_obl

.. autofunction:: u2_u1_obl

.. autofunction:: p2_p1_obl

.. autofunction:: T2_T1_obl

.. autofunction:: p02_p01_obl


sutherland
----------

.. automodule:: cfpylib.gasdyn.sutherland

.. autofunction:: sutherland

.. autofunction:: mu

.. autofunction:: k

