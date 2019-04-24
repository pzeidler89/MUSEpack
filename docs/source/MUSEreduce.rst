.. _musereduce:

MUSEreduce
**********

:class:`MUSEreduce.musereduce` is a easy-to-use python wrapper for the VLT/`MUSE data reduction pipeline`_ and does not replace the core functionalities of the pipeline provided by ESO. In order to function properly we recommend to install the latest pipeline version found under: https://www.eso.org/sci/software/pipelines/muse/.

.. _MUSE data reduction pipeline: https://www.eso.org/sci/software/pipelines/muse/

.. autoclass:: MUSEreduce.musereduce
   :members:

Utility modules
---------------

The following modules are routines to support the :class:`MUSEreduce.musereduce`, needed to execute the various data reduction steps as it is suggested in the `MUSE data reduction pipeline`_. They are documented here for completness.

.. automodule:: MUSEreduce
   :members:
   :private-members:
   :exclude-members: musereduce
   :member-order: bysource

History
-------

.. versionadded:: 0.1.0
   Executes the MUSE reduction pileline in the correct order
.. versionadded:: 0.1.1
   introducing only one calibration file folder per OB
.. versionadded:: 0.1.2
   choosing the illumination file closest to the observation
.. versionadded:: 0.1.3
   selecting the files for the different master file creations
.. versionadded:: 0.1.4
   minor corrections
.. versionadded:: 0.1.5
   looping order changed. each module loops by itself
.. versionadded:: 0.1.6
   always checking if calibration files already exist
.. versionadded:: 0.1.7
   Choosing between ESO calibrations and own calibrations
.. versionadded:: 0.1.8
   User can choose specific exposure time to reduce
.. versionadded:: 0.2.0
   exposures spread via different OBs for one pointing is supported. To do so: run the script normally for each OB including ``muse_scipost``. After all are processed then run ``muse_exp_align`` and ``muse_exp_combine``. AO observations are supported. To reduce multiple OBs in one run, set rootpath to *manual*.
.. versionadded:: 0.2.1
   user can choose the number of CPUs used
.. versionadded:: 0.2.2
   sky subtraction can be modified, so that individual elements can be excluded
.. versionadded:: 0.2.3
   use *json* file as input
.. versionadded:: 0.2.4
   additional parameters added: *skyreject* and *skysubtraction*. The set parameters and modules are shown in the output
.. versionadded:: 0.2.5
   bug fixes
.. versionadded:: 0.3.0
   using the correct *ILLUM* file for the *STD* reduction in the ``muse_sci_basic`` routine, ``muse_sci_basic`` now separated for *STD* and *OBJECT* reduction.
.. versionadded:: 0.3.1
   one can select if the sof file are created automatically or provided by the user
.. versionadded:: 0.4.0
   supports now *pipeline 2.4.2* and the *NFM-AO* added: pipeline_path choosing if darks may be used only reduces *STD* once per OB general use of external *SKY* fields collecting the files for ``muse_exp_combine`` in an independent step ``muse_exp_align`` is an independent step now
.. versionadded:: 0.4.1
   new file names to correct a problem where data gets replaced in the ``muse_scipost`` routine if you reduce the data with and without sky
.. versionadded:: 0.4.2
   one can now change the ignore and fraction parameters in the *json* file
.. versionadded:: 0.4.3
   one can auto remove and rewrite the statics
.. versionadded:: 0.4.4
   changed the sky subtraction keyword the user can give now individual names for the different dither exposures: Does currently not work with multiple OBs or multiple pointings per OB
.. versionadded:: 0.5.0
   rewriting :class:`MUSEreduce.musereduce` to a class and pep-8 style. *DEBUG* keyword added. Wrapper can be executed without running ``esorex`` but needs to be used with already existing reduced data.
.. versionadded:: 0.5.1
   added *skymethod*
.. versionadded:: 0.5.2
   M:class:`MUSEreduce.musereduce` can now handle if the exposures for one pointing are distributed via multiple OBs with multiple exposures in each OB.
