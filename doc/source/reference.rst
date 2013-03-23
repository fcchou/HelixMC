.. _reference:

#################
HelixMC Reference
#################

:Release: |version|
:Date: |today|


.. currentmodule:: helixmc

This reference document details functions, modules, and objects
included in HelixMC, describing what they are and what they do.
For learning how to use HelixMC, see also :ref:`tutorial`.

Constant Random Seed
====================
.. autosummary::
   :toctree: generated/

   constant_seed

Helix Pose
==========
.. autosummary::
   :toctree: generated/
   :template: class.rst

   pose.HelixPose

Random Base-pair Steps Generators
=================================
.. autosummary::
   :toctree: generated/
   :template: class.rst

   random_step.RandomStepBase
   random_step.RandomStepSimple
   random_step.RandomStepAgg

Score Functions
===============
.. autosummary::
   :toctree: generated/
   :template: class.rst

   scorefxn.ScorefxnBase
   scorefxn.ScorefxnTweezers

Utility Functions
=================

Rotation Matrices
-----------------
.. autosummary::
   :toctree: generated/

    util.Rz
    util.Rx
    util.Ry
    util.R_axis

Twist and Writhe
----------------
.. autosummary::
   :toctree: generated/

   util.ribbon_twist
   util.writhe_exact
   util.writhe_fuller

Useful Conversions
------------------
.. autosummary::
   :toctree: generated/

   util.params2coords
   util.coords2params
   util.dr2coord
   util.params2data

Other Functions
---------------
.. autosummary::
   :toctree: generated/

   util.unitarize
   util.MC_acpt_rej
   util.read_seq_from_fasta

Useful Fitting Functions
========================
.. autosummary::
   :toctree: generated/

   fitfxn.wlc_bouchiat
   fitfxn.wlc_bouchiat_impl
   fitfxn.f_wlc_bouchiat_impl
   fitfxn.moroz_3rd
   fitfxn.moroz_1st

