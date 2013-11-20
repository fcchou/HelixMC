Base-pair Step Parameters Database
==================================
In the ``helixmc/data/`` folder, several different bp-step parameter sets are
given. These datasets were all extracted from structures in Protein Data Bank
(PDB, http://www.pdb.org/), with different selection and filtering. The list
below summarizes these data.

:DNA_default:
    B-DNA data from structures with resolution (Rs) <= 2.8 Å,
    excluding protein-binding models.

:DNA_2.8_all:
    A-DNA + B-DNA, Rs <= 2.8 Å, including protein-binding models.

:DNA_2.0_noprot:
    B-DNA, Rs <= 2.0 Å, excluding protein-binding models.

:RNA_default:
    RNA, Rs <= 2.8 Å, excluding protein-binding models.

:RNA_2.8_all:
    RNA, Rs <= 2.8 Å, including protein-binding models.

:RNA_2.0_noprot:
    RNA, Rs <= 2.0 Å, excluding protein-binding models.

:Z-DNA:
    Z-DNA, Rs <= 2.8 Å, including protein-binding models.

:\*unfiltered:
    Unfiltered datasets (no filtering of histogram outliers).

:DNA_gau:
    Single 6D Gaussian built from DNA_default.

:RNA_gau:
    Single 6D Gaussian built from RNA_default.

:DNA_gau_graft:
    Chimera dataset with mean from DNA_gau and covariance from RNA_gau.

:RNA_gau_graft:
    Chimera dataset with mean from RNA_gau and covariance from DNA_gau.

:\*gau_refit
    Manually refitted datasets to match experimental measurements.

:\*_2.8_all_\?bp:
    Multi-bp datasets derived from the `2.8_all` pdb lists.

Note that Gaussian dataset (`*gau*.npy`) must be loaded with
`-gaussian_params` tag in `helixmc-run` command line (instead of `-params`).
Also Gaussian dataset does not support sequence specific simulations.

The corresponding lists of PDB models being used are given in the
``helixmc/data/pdb_list/`` folder.
