.. MUSEpack documentation master file, created by
   sphinx-quickstart on Fri Mar  1 14:10:22 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MUSEpack - a tool to analyze IFS datasets
=========================================

Introduction
------------

MUSEpack is a python package written to support the data analyzes from Integral Field Units, specifically tailored to use datasets of the Multi Unit Spectroscopic Explorer (`MUSE`_) mounted at UT4 of the VLT.

The main purpose of MUSEpack is the measurement of stellar and gas radial velocities to an accuracy of :math:`1-3\,\rm{km/s}` without the need for spectral template libraries using :class:`radial_velocities.RV_spectrum`. By using strong stellar absorption lines or gas emission lines and a high-resolution spectral line library such as the `NIST Atomic Spectra Database`_ it is possible to create template, which are fitted to the spectra using a Monte Carlo method. For a detailed description of the code we refer to `Zeidler et al. 2019`_.

This package includes the class :class:`MUSEreduce.musereduce`, which is a easy-to-use python wrapper for the `MUSE data reduction pipeline`_ provided by ESO including a modification to ``muse_create_sky`` which prevents over subtraction of emission lines when reducing data in regions with a high, variable background such as Galactic HII regions.


.. toctree::
   :maxdepth: 2
   :numbered:
   :caption: Contents:

   installation
   MUSEreduce
   radialvelocities
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: Disclaimer:

   disclaimer

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`


.. _MUSE: https://www.eso.org/sci/facilities/develop/instruments/muse.html
.. _MUSE data reduction pipeline: https://www.eso.org/sci/software/pipelines/muse/
.. _NIST Atomic Spectra Database: https://www.nist.gov/pml/atomic-spectra-database
.. _Zeidler et al. 2019: www.xyz.com


References
----------

* Bacon, R., Accardo, M., Adjali, L., et al. 2010, Proc. SPIE, 7735, 773508
* Kamann, S., Husser, T.-O., Brinchmann, J., et al. 2016, A&A, 588, A149
