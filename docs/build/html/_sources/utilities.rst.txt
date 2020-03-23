.. _utility module:

Utility module
**************

This modules contains a variety of functions to support *MUSEpack*. Many of these functions may be useful for other purposes.

.. todo::

   The documentation of all of the utility modules will follow soon

.. automodule:: utils
   :members:


.. _Zeidler et al. 2019: https://ui.adsabs.harvard.edu/abs/2019AJ....158..201Z/abstract

.. _Cappellari and Emsellem 2004: https://ui.adsabs.harvard.edu/abs/2004PASP..116..138C/abstract


History
-------

.. versionadded:: 0.1.0
   module created

.. versionadded:: 0.1.1
   moved to pep-8

.. versionadded:: 0.1.2
   now handles absorption and emission lines emission not tested yet, though

.. versionadded:: 0.1.3
   The :mod:`util.Line_clipping` was adjusted in how the two outliers are clipped before the MAD is calculated. It now keeps the :math:`N-2`-lines
   that have the smaller deviation from each other.

.. versionadded:: 0.1.4
   adding `rv_sys`, to compensate for larger systematic RV shifts or redshifts for the line and RV fitter. The :mod:`util.lambda_rv_shift` was introduced