MUSEpack v1.0.2 (legacy)
========================

.. warning::
   **This version of MUSEpack v1.0.2 will be kept as a legacy version that specifically supports ESO's MUSE data reduction v2.4.2. This particular version can be    downloaded from https://github.com/pzeidler89/MUSEpack/tree/legacy_MUSEpipeline_v2.4.2**

Introduction
------------
*MUSEpack* is a python package written to support the data analyzes from Integral Field Units, specifically tailored to use datasets of the Multi Unit Spectroscopic Explorer (`MUSE`_) mounted at UT4 of the VLT (`Bacon et al. 2010`_).

The main purpose of MUSEpack is to measure stellar and gas radial velocities to an accuracy of :math:`1-2\,{\rm kms}\,{\rm s}^{-1}` without the need for spectral template libraries using :class:`radial_velocities.RV_spectrum`. With strong stellar absorption lines or gas emission lines and a high-resolution spectral line library such as the `NIST Atomic Spectra Database`_ it is possible to create templates using the observed spectra, which can be cross correlated with the spectra using a Monte Carlo method. For a detailed description and the citation of the code we refer to `Zeidler et al. 2019`_.

On this website we introduce the individual python classes and modules of *MUSEpack*. In :ref:`examples` we present detailed instructions on how to use the main modules and classes of *MUSEpack*. For any suggestions and bug reports please create a pull request in GitHub https://github.com/pzeidler89/MUSEpack or send an email to peterzeidler89@gmail.com.


.. image:: https://zenodo.org/badge/144210611.svg
   :target: https://zenodo.org/badge/latestdoi/144210611



.. toctree::
   :maxdepth: 2
   :numbered:
   :caption: Contents:

   installation
   MUSEreduce
   radialvelocities
   cubes
   utilities
   examples

.. toctree::
   :maxdepth: 2
   :caption: Credits:

   credits

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
.. _Zeidler et al. 2019: https://ui.adsabs.harvard.edu/abs/2019AJ....158..201Z/abstract
.. _Bacon et al. 2010: https://ui.adsabs.harvard.edu/abs/2010SPIE.7735E..08B/abstract

References
----------

* Zeidler, P., Nota, A., Sabbi, E., Luljak, P., McLeod, A. F., Grebel, E. K., Pasquali, A., Tosi, M. 2019, AJ, 158, 201
* Bacon, R., Accardo, M., Adjali, L., et al. 2010, Proc. SPIE, 7735, 773508