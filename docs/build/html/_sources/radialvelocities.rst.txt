.. _radial velocities:

Spectral line fitting and radial velocity measurements
******************************************************

The main purpose of MUSEreduce is the spectral line fitting in order to create template spectra from observations, which are used to measure the radial velocities of stars and gas, especially in the absence of working spectral template libraries, e.g., for pre-main-sequence stars in the optical.

The :class:`radial_velocities.RV_spectrum` class and the spectral fitting module :mod:`line_fitter` are the heart of this package using `pyspeckit`_ as basic spectral fitting routine and `ppxf`_ for the spectral cross-correlation. Together with a Monte Carlo bootstrap technique, it is possible to measure radial velocities with an accuracy of :math:`1-3\,\rm{km}\,\rm{s}^{-1}` in the absence of any template library. For a detailed description on how this technique works, we refer to `Zeidler et al. 2019`_.

.. _pyspeckit: https://pyspeckit.readthedocs.io/en/latest/index.html
.. _ppxf: http://www-astro.physics.ox.ac.uk/~mxc/software/#ppxf


Spectral line fitter
--------------------

This module is the work horse for fitting spectral lines from an input catalog to a given 1D spectrum. The line fitting itself is performed using the module `pyspeckit`_.

`pyspeckit`_ uses an iterative method to fit the spectral lines and the continuum simultaneously. :mod:`line_fitter` automatically adjusts the input parameters to `pyspeckit`_ in the manner to optimize the fitting (*continuum order* and *wavelegth limits*)

The :mod:`line_fitter` determines whether the fit fullfils the set parameters and was succesfull. The :mod:`line_fitter` also allows to fit blends by fixing a maximum ratio between the primary line (the one the user is interested in) and the secondary line (blend). This ensures that the blend does not become the dominant line.

In order to accomodate hyper-velocity stars an option is given that the wavelength limits are automatically adjust for each iteration :math:`n` based on the solution :math:`n-1` of the primary line taking into account :math:`\Delta \lambda / \lambda`.

If the spectra with high radial velocities, or if the object is red shifted an initial guess of the radial velocity needs to be provided. This may be done by the *Kwargs* `rv_sys` of the :class:`radial_velocities.RV_spectrum`. This is important if the shift is much larger than the wavelength limits, where the automatic adjustment fails. Additionally, the more accurate this `rv_sys` is provided the faster the fit converges.

.. _pyspeckit: https://pyspeckit.readthedocs.io/en/latest/index.html

.. automodule:: line_fitter
   :members:

.. warning::

   The automatic handling of absorption and emission lines implemented in vers. 0.1.2 has not been tested yet.

Radial Velocities
-----------------

This is the main class for measuring the radial velocities. A basic example is descibed in :ref:`examples`. A detailed description can be found in Zeidler et al. 2019.

Throughout the document the **primary line** is the spectral line of interest for which the line fitting will be executed. **Secondary lines** are spectral lines in blends.

A detailed demonstration on how to use the radial velocity fitter is provided in :ref:`examples` including template files.

The radial velocity fitter
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: radial_velocities.RV_spectrum
   :members:

.. _Monte Carlo:

Monte Carlo radial velocity determination
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is the module that performs the Monte Carlo boot strapping to measure radial velocities. It is automatically called by the :class:`radial_velocities.RV_spectrum.rv_fit` attribute of :class:`radial_velocities.RV_spectrum` class and there is no need to repeat this step manually. Nevertheless, we show this part of the code since it may be useful for other applications to measure RVs.

.. automodule:: ppxf_MC
   :members:


.. _Zeidler et al. 2019: www.xyz.com

.. _Cappellari and Emsellem 2004: https://ui.adsabs.harvard.edu/abs/2004PASP..116..138C/abstract

History
-------
.. versionadded:: 0.1.0
   working radial velocity fitting package

.. versionadded:: 0.1.0
   no RV calculation for lines that are below the significance level
.. versionadded:: 0.1.0
   instrument dispersion added as keyword
.. versionadded:: 0.1.1
   moved to pep-8
.. versionadded:: 0.1.2
   now handles absorption and emission lines. **Not tested yet, though**
.. versionadded:: 0.1.3
   During the boot strap procedure, the initial guesses of the velocities are varied between certain limits, to ensure more stability to the
   fit in the case of an "unlucky choice of initial parameters. Added Kwarg to :mod:`RV_spectrum.rv_fit` is `RV_guess_var`, which is by default 0. 
   It describes the min/max variation of the initial RV guess for each fit.
.. versionadded:: 1.0
   The release version as originally published in `Zeidler et al. 2019`_.
.. versionadded:: 1.1
   Introducing the *Kwargs* `rv_sys`, which should be used for large initial RV offsets from rest frame or for redshifted spectra. `rv_sys` should provide an estimate of that shift.