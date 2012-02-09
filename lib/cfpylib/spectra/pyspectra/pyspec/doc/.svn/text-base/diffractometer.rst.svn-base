Diffractometer
==============

This class provides the transformation from the setting angles of a six circle diffractometer into reciprocal space. There for we use mainly four frames in reciprocal Space: :math:`\theta` -frame, :math:`\phi` -frame, :math:`U` -frame and :math:`(H,K,L)` -frame.

Lab Frame
---------

We define the lab frame :math:`\Sigma` like in the paper about the four circle diffraction setup (Busing1967), i.e. the :math:`Z_4` -axis up, the :math:`Y_4` -axis along the X-ray beam and the :math:`X_4` -axis forms a right handed orthonganl system with the others. The figure shows the setup with all six angles and a diffrent lab frame from the paper Lohmeier1993 (:math:`X_6=Z_4`, :math:`Y_6=Y_4`, :math:`Z_6=-X_4`).

.. math::

   \bf{a}_4 = V_{46} \bf{a}_6 = \begin{pmatrix} 0 & 0 & -1 \\
                                                0 & 1 &  0 \\
		 		                1 & 0 &  0 \end{pmatrix} \bf{a}_6

.. image:: figures/Lohmeier1993_Fig1.*

Rotations and Setting Angles
----------------------------

The rotation by :math:`\mu` along the :math:`+Z_4`\ -axis leads to :math:`\Sigma'` followed by the angles familiar from the four circle setup, i.e.
rotation by :math:`\theta` along the :math:`+X_4'`\ -axis leads to :math:`\Sigma''` 
(:math:`\theta`\ -frame, :math:`\Sigma^\theta`),
rotation by :math:`\chi`   along the :math:`+Y_4''`\ -axis leads to :math:`\Sigma'''`, and
rotation by :math:`\phi`   along the :math:`+X_4'''`\ -axis leads to :math:`\Sigma''''`
(:math:`\phi`\ -frame, :math:`\Sigma^\phi`)
The detector rotations start in :math:`\Sigma'` with 
rotation by :math:`\delta` along the :math:`+X_4'` -axis into :math:`\Sigma^*`, and
rotation by :math:`\gamma` along the :math:`+Z_4^*` -axis into :math:`\Sigma^{**}`.
Note that in the origanal paper :math:`\mu` is called :math:`\alpha`.
After this introduction we will drop most times the indices at the axes and refer to the four circle lab frame.

.. table:: Rotations and frames for six circle geometry

   ==============  ===============  ===============  ==========================================
   Angle           Axis\ :sub:`4`   Axis\ :sub:`6`   New Frame
   ==============  ===============  ===============  ==========================================
   :math:`\mu`     :math:`+Z_4`     :math:`+X_6`     :math:`\Sigma'`
   :math:`\theta`  :math:`+X_4'`    :math:`-Z_6'`    :math:`\Sigma''`   (:math:`\Sigma^\theta`)
   :math:`\chi`    :math:`+Y_4''`   :math:`+Y_6''`   :math:`\Sigma'''`
   :math:`\phi`    :math:`+X_4'''`  :math:`-Z_6'''`  :math:`\Sigma''''` (:math:`\Sigma^\phi`)
   :math:`\delta`  :math:`+X_4'`    :math:`-Z_6'`    :math:`\Sigma^*`
   :math:`\gamma`  :math:`+Z_4^*`   :math:`+X_6^*`   :math:`\Sigma^{**}`
   ==============  ===============  ===============  ==========================================

Rotation matrcies
-----------------

Rotation of a vector :math:`\bf{a}` along the :math:`Z`\ -axis by the angle :math:`\alpha` is discribed by the rotaion matrix :math:`R_Z(\alpha)`

.. math::

   R_Z(\alpha) &= \begin{pmatrix} \cos(\alpha) & -\sin(\alpha) & 0 \\
                                  \sin(\alpha) &  \cos(\alpha) & 0 \\
		 		  0            &  0            & 1 \end{pmatrix} \\
   \bf{a}' &= R_Z(\alpha) \bf{a}	 

The back transformation is the transposed :math:`R_Z(\alpha)^\mathrm{T}` or :math:`R_Z(-\alpha)`, the back rotation. If the coordinate system is rotated by the angle :math:`\alpha`, :math:`R_Z(-\alpha)` is the transformation.
For completness the other two rotation matrices along the principle axes.

.. math::

   R_X(\alpha) &= \begin{pmatrix} 1 & -\sin(\alpha) &  0            \\
                                  0 &  \cos(\alpha) & -\sin(\alpha) \\
		 	          0 &  \sin(\alpha) &  \cos(\alpha) \end{pmatrix} \\
   R_Y(\alpha) &= \begin{pmatrix}  \cos(\alpha) &  0 & \sin(\alpha) \\
                                   0            &  1 & 0            \\
		         	  -\sin(\alpha) &  0 & \cos(\alpha) \end{pmatrix}


Calculation of :math:`Q`
------------------------

The initial and final (scattered) wavevectors are :math:`\bf{k}_\mathrm{i}`, :math:`\bf{k}_\mathrm{f}` and have the same length :math:`k` (elastic scattering).
For all setting angels set to 0, :math:`\bf{k}_\mathrm{i}` and :math:`\bf{k}_\mathrm{f}` lay along the :math:`Y`\ -axis.

.. math::
  
   \bf{k}_i(0) = \bf{k}_f(0) = \begin{pmatrix} 0 \\ k \\ 0 \end{pmatrix}

To get :math:`\bf{k}_\mathrm{i}''`, i.e. in :math:`\Sigma_\theta`, we have to apply the rotation matrix along :math:`Z` by the angle :math:`-\mu` and the rotation matrix along :math:`X'` by the angle :math:`-\theta`, because :math:`\bf{k}_\mathrm{i} = k_i(0)` for all angles. To get :math:`\bf{k}_\mathrm{f}''` we have to considered that :math:`\bf{k}_\mathrm{f}^** = k_f(0)` for all angles, becuase the measured scattered X-rays hit directly the detector. In a similar sence we need rotations by :math:`+\gamma` along :math:`Z^*` and by :math:`-\theta + \delta` along :math:`X'`.

.. math::
  
   \bf{k}_i'' &= R_X(-\theta) R_Z(-\mu) \bf{k}_i(0) \\
   \bf{k}_f'' &= R_X(-\theta + \delta) R_Z(\gamma) \bf{k}_f(0)

Now we get :math:`\bf{Q}^\theta` in :math:`\Sigma^\theta`.

.. math::

   \bf{Q}^\theta = \bf{Q}'' = \bf{k}_f'' - \bf{k}_i''

If we take the rotations by :math:`\chi` and :math:`\phi` into account we come to :math:`\Sigma^\phi`.

.. math::

   \bf{Q}^\phi = R_Z(-\phi) R_Y(-\chi) \bf{Q}^\theta

The last step is transforming from :math:`\Sigma^\phi_4` to :math:`\Sigma^\phi_6` and applying the inverse orientation matrix :math:`UB^{-1}` from e.g. Spec SIXC to get the :math:`HKL`\ -values.

.. math::

   \bf{Q}^{HKL} = UB^{-1} V_{46}^\mathrm{T} \bf{Q}^\phi

Diffractometer Class
--------------------

.. automodule:: pyspec.diffractometer
   :members:
