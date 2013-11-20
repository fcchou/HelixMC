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

Code Organization
=================
HelixMC is coded in Python in an object-oriented fashion, allowing easy usage
and extension. Here we briefly summarizes the organization of the code.

First, :ref:`HelixPose <HelixPose>` object stores all information of the
current helix conformation. Various properties, e.g. coordinates, twist, writhe,
etc. can be directly accessed from the object. The object also contains
functions that allow one to update the conformation and to plot the helix.

Second, :ref:`RandomStep <RandomStep>` are used to generate random samples of
base-pair step parameters. `RandomStepSimple` takes in a list of database
parameter sets, and can randomly emit one of input parameter set, or construct
a multivariate Gaussian and emit samples from the distribution, upon the choice
of the user. `RandomStepAgg` aggregates multiple `RandomStep` objects into one,
allow easy handling of sequence-dependent sampling (by aggregating several
`RandomStepSimple` for different sequences). The user can also create their
own `RandomStep` objects by inheriting from `RandomStepBase`.

Third, :ref:`Score <Score>` objects take a `HelixPose` and scores it.
The scoring can then be used to decide whether a Monte-Carlo move should be
accepted. Currently we have 3 simple score terms `ScoreExt`,
`ScoreTorsionTrap` and `ScoreXyTrap`, which scores a HelixPose under
Z-extension, torsional trap and xy horizontal trap respectively. These 3 score
terms are summarized into `ScoreTweezers`, which is a sub-class of `ScoreAgg`,
an aggregator class that combines multiple score functions. `ScoreTweezers` is
the workhorse score functions currently in used. The user can also define
their own scoring by inheriting from `ScoreBase`.

Last, the :ref:`util <util>` module contains useful functions for evaluating
twist and writhe, for conversion between bp-step parameters and cartesian
translation and rotation operation, and so on. The :ref:`fitfxn <fitfxn>`
module is a standalone module that contains a few widely used analytical
fitting functions based on the elastic rod model.

Constant Random Seed
====================
.. autosummary::
   :toctree: generated/

   constant_seed

.. _HelixPose:

Helix Pose
==========
.. autosummary::
   :toctree: generated/
   :template: class.rst

   pose.HelixPose

.. _RandomStep:

Random Base-pair Steps Generators
=================================
.. autosummary::
   :toctree: generated/
   :template: class.rst

   random_step.RandomStepBase
   random_step.RandomStepSimple
   random_step.RandomStepAgg

.. _Score:

Score Functions
===============
.. autosummary::
   :toctree: generated/
   :template: class.rst

   score.ScoreBase
   score.ScoreExt
   score.ScoreTorsionTrap
   score.ScoreXyTrap
   score.ScoreAgg
   score.ScoreTweezers

.. _util:

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
   util.dr2coords
   util.coords2dr
   util.params2data
   util.data2params
   util.params_join
   util.frames2params_3dna

Other Functions
---------------
.. autosummary::
   :toctree: generated/

   util.unitarize
   util.MC_acpt_rej
   util.read_seq_from_fasta

.. _fitfxn:

Useful Fitting Functions
========================
.. autosummary::
   :toctree: generated/

   fitfxn.wlc_odijk
   fitfxn.wlc_bouchiat
   fitfxn.wlc_bouchiat_impl
   fitfxn.f_wlc_bouchiat_impl
   fitfxn.moroz_3rd
   fitfxn.moroz_1st

