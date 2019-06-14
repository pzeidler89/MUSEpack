.. _radial velocities:

Spectral line fitting and radial velocity measurements
******************************************************

The main purpose of MUSEreduce is the spectral line fitting in order to create template spectra from observations, which are used to measure the radial velocities of stars and gas, especially in the absence of working spectral template libraries, e.e., for pre-main-sequence stars.

The :class:`radial_velocities.RV_spectrum` and the spectral fitting module :mod:`line_fitter` are the heart of this package. By using `pyspeckit`_ as basic spectral fiutting routine and `ppxf`_ for the spectral cross-correlation, including a Monte Carlo bootstrap technique, it is possible to measure radial velocities with an accuracy of :math:`1-3\,\rm{km/s}` in the absence of any template library. For a detailed description we refer to Zeidler et al. 2019.

.. _pyspeckit: https://pyspeckit.readthedocs.io/en/latest/index.html
.. _ppxf: http://www-astro.physics.ox.ac.uk/~mxc/software/#ppxf


Radial Velocities
-----------------

This is the main class for measuring the radial velocities. A basic example is descibed in :ref:`examples`. A detailed description can be found in Zeidler et al. 2019.

Throughout the document the **primary line** is the spectral line of interest for which the line fitting will be executed. **Secondary lines** are spectral lines in blends.

.. warning::

   The automatic handling of absorption and emission lines since vers. 0.1.2 has not been tested yet.

.. autoclass:: radial_velocities.RV_spectrum
   :members:

Monte Carlo radial velocity determination
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is the module that does the Monte Carlo boot strapping to measure RV to desired accuracy. It is automatically called by the ``rv_fit`` attribute of :class:`radial_velocities.RV_spectrum` and there is no need to repeat this step manually. Nevertheless, we show this part of the code since it may be useful for other applications to measure RVs.

.. automodule:: ppxf_MC
   :members:


Spectral line fitter
--------------------

This module is the work horse for fitting spectral lines from an input catalog to a given 1D spectrum. The line fitting itself is performed using the module `pyspeckit`_.

`pyspeckit`_ uses an iterative method to fit the spectral lines and the continuum simultaneously. :mod:`line_fitter` automatically adjusts the input parameters to `pyspeckit`_ in the manner to optimize the fitting (*continuum order* and *wavelegth limits*)

The :mod:`line_fitter` determines whether the fit fullfils the set parameters and was, therefore, succesfull and can be further used to do the velocity fitting. :mod:`line_fitter` also allows to fit blends by fixing a maximum ratio between the primary line (the one the user is interested in) and the secondary line (blend). This ensures that the blend does not become the dominant line.

In order to accomodate hyper-velocity stars an option is given that the wavelength limits are automatically adjust for each iteration based on the solution of the primary line.

.. _pyspeckit: https://pyspeckit.readthedocs.io/en/latest/index.html

.. automodule:: line_fitter
   :members:

Cross correlating the spectra
-----------------------------

.. todo::

   The documentation will follow soon

Utility modules
---------------

.. todo::

   The documentation of utility modules will follow soon




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