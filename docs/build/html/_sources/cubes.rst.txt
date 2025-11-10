.. _cubes:

Cubes
*****

:mod:`cubes` is a module that contains a collection of support functions for the analyses of data cubes, specifically MUSE data cubes.

.. todo::

   A more detailed description will follow soon.

.. automodule:: cubes
   :members:

History
-------

.. versionadded:: 0.1.0
   The functions :func:`cubes.wcs_cor`, :func:`cubes.pampelmuse_cat`, :func:`cubes.linemaps`, and :func:`cubes.mosaics` were added.

.. versionadded:: 0.1.0
   :func:`cubes.wcs_cor` also can be used with *the OFFSET_LIST.fits* file produced by ESORex

.. versionadded:: 0.2.0
    :func:`cubes.wcs_cor` can now also be used with any zeropoints of any filter. This is especially important if we match to Gaia DR3 directly

.. versionadded:: 0.2.0
    minor bugfixes