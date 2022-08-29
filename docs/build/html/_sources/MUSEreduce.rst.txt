.. _musereduce:

MUSEreduce
**********

:class:`MUSEreduce.musereduce` is an easy-to-use python Class used as wrapper for the VLT/`MUSE data reduction pipeline`_ and does not replace the core functionalities of the pipeline provided by ESO. In order to function properly with this version of :class:`MUSEreduce.musereduce` we recommend to install the pipeline version v2.8.7 found under: https://www.eso.org/sci/software/pipelines/muse/. This version of :class:`MUSEreduce.musereduce` has been tested with version v2.8.1 of the data reduction pipeline.

.. _MUSE data reduction pipeline: https://www.eso.org/sci/software/pipelines/muse/

To run :class:`MUSEreduce.musereduce` all of the raw data must be stored in a folder named `user_path/raw/OB_ID`, where `OB_ID` must be a unique name for each individual OB. It is not mandatory (but recommended) to use the same nomenclature as as in the fits header. The script has to point to `user_path` by setting the keyword ``rootpath`` of the :obj:`json` config file to the parent directory.

.. autoclass:: MUSEreduce.musereduce
   :members:

The basic folder structure
--------------------------

If ``auto_sort_data`` = :obj:`True` in the config file and :class:`MUSEreduce.musereduce` is executed, the basic folder structure is created in the ``rootpath``. As example, we consider that three OBs were observed: OB1a, OB1b, and OB2. OB1a and OB1b are two dither positions (with rotation angles of 5 and 95 deg and an exposure time of 2800s) of the same pointing OB1, while OB2 consists of three dither postions observed in one OB (with rotation angles of 10, 100, and 190 deg and an exposure time of 220s). The folder structure looks the following:

.. image:: images/folder_struct.png

Each *master* OB (`OB1`, `OB2`) has its own folder in the `reduced` folder. Pointing `OB1` consists of the two OBs `OB1a` and `OB1b`. Each OB has the following folders:

 * Each individual exposure (e.g., `091315-402023_2800_005_00`) following the structure: **RADEC_EXPTIME_ROTANGLE_COUNTER**. **RADEC** are the coordinates in sexagesimal, **EXPTIME** is the exposure time in seconds and the **ROTANGLE** is the rotation angle in degrees. The **COUNTER** is needed if there are two exposures with the same configuration.
 * Each individual pointing (e.g., `091315-402023_2800`) following the structure: **RADEC_EXPTIME**.
 * The `calibrations` contain `DARK`, `TWILIGHT`, and `SCIENCE`, in case the calibration files are created from the raw calibration files provided by ESO Therefore, the calibration steps must be executed.
 * The `ESO_calibrations` is the folder, into which the reduced calibration files delivered by ESO (if available) are copied.
 * The `static_calibration_files`, is the folder, into which the statics (part of the `MUSE data reduction pipeline`_ installation) are copied.
 * The folder `std` will contain the reduced data of the standard star.

.. note::
   The calibration files and the standard star are the same per individual OB and are only needed once.

The config file
---------------

The :obj:`json` configuration file is needed to run :class:`MUSEreduce.musereduce`. This file contains all of the data reduction setup. The :obj:`json` configuration file can be found in the main directly of MUSEpack and can be downloaded :download:`here <../../MUSEpack/config.json>`.

.. literalinclude:: ../../MUSEpack/config.json
  :language: JSON

The file :download:`config file <../../MUSEpack/config.json>` is structured the following. If keywords are directly controlling toggles of the `MUSE data reduction pipeline`_ their naming is identical.


.. json:object:: config

   :property global: global parameters affecting the data reduction process
   :proptype global: :json:object:`global`
   
   :property calibration: parameters that affect the calibration steps of the `MUSE data reduction pipeline`_
   :proptype calibration: :json:object:`calibration`
   
   :property sci_basic: the preprocessing of the science exposures
   :proptype sci_basic: :json:object:`sci_basic`
   
   :property std_flux: the flux calibration
   :proptype std_flux: :json:object:`std_flux`
   
   :property sky: the sky subtraction
   :proptype sky: :json:object:`sky`
   
   :property sci_post: the postprocessing step of the data reduction
   :proptype sci_post: :json:object:`sci_post`
   
   :property dither_collect: collecting the individual data cubes and pixel tables to combine them to one pointing
   :proptype dither_collect: :json:object:`dither_collect`
   
   :property exp_align: aligning the individual data cubes and pixel tables to one common WCS
   :proptype exp_align: :json:object:`exp_align`
   
   :property exp_combine: combining the individual data cubes and pixel tables to one final science product
   :proptype exp_combine: :json:object:`exp_combine`

.. json:object:: global

   :property pipeline_path: The absolut path to the `MUSE data reduction pipeline`_ installation folder.
   :proptype pipeline_path: string

   :property mode: The observation mode the data was obtained with.
   :proptype mode: string
   :options mode: WFM-NOAO, WFM-AO, NFM-AO

   :property withrvcorr: bariocentric correction. Needs to be turned off, if one wants run an own wavelength calibration
   :proptype withrvcorr: bool
   :options withrvcorr: true, false, default='true'

   :property auto_sort_data: The raw data is sorted and the calibration files are assigned based on their header information. If :obj:`True`,
                             the file lists `ID_DAR.list`, `ID_SCI.list`, and `ID_TWI.list` are created. If these have to be altered manually (e.g., using
                             different calibration files), we recommend to run it first with ``auto_sort_data`` = :obj:`True`, then make the changes accordingly
                             and from this point on set ``auto_sort_data`` = :obj:`False`.
   :proptype auto_sort_data: bool
   :options auto_sort_data: false, true, default=true

   :property using_specific_exposure_time: The user can choose to only reduce a specific exposure time, if the same OB contains multiple exposures with
                                            different exposure times (e.g., long and short exposures)
   :proptype using_specific_exposure_time: float

   :property dither_multiple_OBs: Each OB is normally limited to a total exposure time of one hour. Therefore, one pointing may be distributed via multiple OBs.
                                 If ``dither_multiple_OBs`` = :obj:`True` it is possible to dither exposures from multiple OBs.
                                 In this case one must provide an ``OB_list``.
   :proptype dither_multiple_OBs: bool
   :options dither_multiple_OBs: false, true, default=false

   :property n_CPU: The number of CPUs used to reduce the data. If set to -1 all available cores are used.
   :proptype n_CPU: int
   :options n_CPU: default=-1

   :property rootpath: The absolut path in which the `raw` folder is located and in which the `processed` folder will be created.
   :proptype rootpath: string

   :property OB: The name of the OB that shall be reduce. It must identical to the OB_ID given in the `raw` folder.
   :proptype OB: string

   :property OB_list: If ``dither_multiple_OBs`` = :obj:`True` one must give a :obj:`string` list of OB_IDs, which will be dithered in the end.
                     The calibration and data reduction runs on each individual OB to take into account different calibration files.
   :proptype OB_list: list
   :options OB_list: default=[]

.. json:object:: calibration

   :property execute: Must be set :obj:`True` if this step should be execute
   :proptype execute: bool
   :options execute: false, true, default=true

   :property using_ESO_Calibration: If the user wants to use the ESO calibration files (recommended if available) instead of running the calibration themselves (BIAS,                                           DARK, FLAT, WAVECAL, LSF, TWILIGHT) it must be set to :obj:`True`
   :proptype using_ESO_Calibration: bool
   :options using_ESO_Calibration: false, true, default=true

   :property renew_statics: Set to :obj:`True` if the static calibration files should be copied again from the `MUSE data reduction pipeline`_ folder. This is necessary if
                             the installed version of the pipline changes and one wants to obtain the latest static calibrations
   :proptype renew_statics: bool
   :options renew_statics: false, true, default=false

   :property dark: Set to :obj:`True` if one wants to also use the DARK files from the calibration.
   :proptype dark: bool
   :options dark: false, true, default=false

   :property esorex_kwargs_bias: A string of additional parameters for the `BIAS` recipe from the MUSE Pipeline User Manual that are not included in the :obj:`json` configuration file otherwise can be provided
   :proptype esorex_kwargs_bias: string
   :options esorex_kwargs_bias: default='none'

   :property esorex_kwargs_dark: A string of additional parameters for the `DARK` recipe from the MUSE Pipeline User Manual that are not included in the :obj:`json` configuration file otherwise can be provided
   :proptype esorex_kwargs_dark: string
   :options esorex_kwargs_dark: default='none'

   :property esorex_kwargs_flat: A string of additional parameters for the `FLAT` recipe from the MUSE Pipeline User Manual that are not included in the :obj:`json` configuration file otherwise can be provided
   :proptype esorex_kwargs_flat: string
   :options esorex_kwargs_flat: default='none'

   :property esorex_kwargs_wavecal: A string of additional parameters for the `WAVECAL` recipe from the MUSE Pipeline User Manual that are not included in the :obj:`json` configuration file otherwise can be provided
   :proptype esorex_kwargs_wavecal: string
   :options esorex_kwargs_wavecal: default='none'

   :property esorex_kwargs_lsf: A string of additional parameters for the `LSF` recipe from the MUSE Pipeline User Manual that are not included in the :obj:`json` configuration file otherwise can be provided
   :proptype esorex_kwargs_lsf: string
   :options esorex_kwargs_lsf: default='none'

   :property esorex_kwargs_twilight: A string of additional parameters for the `TWILIGHT` recipe from the MUSE Pipeline User Manual that are not included in the :obj:`json` configuration file otherwise can be provided
   :proptype esorex_kwargs_twilight: string
   :options esorex_kwargs_twilight: default='none'

   :property create_sof: Must be set :obj:`True` if one wants to create a new `.sof` file. In case the user wants to change or create the `.sof` file manually,
                        than it must be set to :obj:`False` since it will be overwritten otherwise
   :proptype create_sof: bool
   :options create_sof: false, true, default=true

.. json:object:: sci_basic

   :property execute: Must be set :obj:`True` if this step should be executed
   :proptype execute: bool
   :options execute: false, true, default=true

   :property execute_std: Must be set :obj:`True` if the STD star should be reduced
   :proptype execute_std: bool
   :options execute_std: false, true, default=true

   :property sky_reject: The sigma clipping parameters for the Gaussian fit to each sky emission line: `high sigma clipping limit` (:obj:`float`),
                         `low sigma clipping limit` (:obj:`float`), `number of iterations` (:obj:`int`). For a more detailed description we refer to the
                         `MUSE data reduction pipeline`_ manual.
   :proptype sky_reject: string
   :options sky_reject: default='15.,15.,1'

   :property skylines: The sky lines used to calibrate the wavelength offset in each IFU. the default are two lines but more strong skylines may be provided. For a
                       more detailed description we refer to the `MUSE data reduction pipeline`_ manual.
   :proptype skylines: string
   :options skylines: default='5577.339,6300.304'

   :property esorex_kwargs: A string of additional parameters for the `SCIBASIC` recipe from the MUSE Pipeline User Manual that are not included in the :obj:`json` configuration file otherwise can be provided
   :proptype esorex_kwargs: string
   :options esorex_kwargs: default='none'

   :property create_sof: Must be set :obj:`True` if one wants to create a new `.sof` file. In case the user wants to change or create the `.sof` file manually,
                        than it must be set to :obj:`False` since it will be overwritten otherwise
   :proptype create_sof: bool
   :options create_sof: false, true, default=true

.. json:object:: std_flux

   :property execute: Must be set :obj:`True` if this step should be executed
   :proptype execute: bool
   :options execute: false, true, default=true

   :property esorex_kwargs: A string of additional parameters for the `STANDARD` recipe from the MUSE Pipeline User Manual that are not included in the :obj:`json` configuration file otherwise can be provided
   :proptype esorex_kwargs: string
   :options esorex_kwargs: default='none'

   :property create_sof: Must be set :obj:`True` if one wants to create a new `.sof` file. In case the user wants to change or create the `.sof` file manually,
                        than it must be set to :obj:`False` since it will be overwritten otherwise
   :proptype create_sof: bool
   :options create_sof: false, true, default=true

.. json:object:: sky

   :property execute: Must be set :obj:`True` if this step should be executed
   :proptype execute: bool
   :options execute: false, true, default=true

   :property modified: If set :obj:`True` the modified sky subtraction will be executed. This method will prevent the over subtraction of emission lines that are both emitted from the Earth's atmosphere (tellurics) and the target (e.g., HII regions). For a detailed description of this method we refer to `Zeidler et al. 2019`_. If ``modified`` = :obj:`True` the ``method`` should be `subtract_model`. Other methods will **not** lead to an error but they may lead to wrong results.
   :proptype modified: bool
   :options modified: false, true, default=false

   :property sky_field: Determines if a sky observation (if available) or the science observation itself is used to determine the background
                        contamination (tellurics). If set to `auto` the pipeline will check if there are sky observations available and use the closest one in time to
                        the science exposures. If the ``skyfield``  = `object`, the science exposure will be used.
   :proptype sky_field: string
   :options sky_field: 'auto', 'object', default='auto'

   :property fraction: The fraction of the image used to be considered as sky. Must be between 0. and 1. For more details we refer to the
                         `MUSE data reduction pipeline`_ handbook.
   :proptype fraction: float
   :options fraction: [0.,1.[, default=0.05

   :property ignore: The fraction of the image ignored for sky consideration. Must be between 0. and 1. For more details we refer to the
                         `MUSE data reduction pipeline`_ handbook.
   :proptype ignore: float
   :options ignore: [0.,1.[, default=0.05

   :property method: The method how the determined sky spectrum is fitted to each spaxel to subtract the sky background. Should be set to 'subtract-model' in case
                     of ``modified`` = :obj:`True`. For more details we refer to the `MUSE data reduction pipeline`_ handbook.
   :proptype method: string
   :options method: 'model', 'subtract-model', 'simple', 'none', default='model'

   :property esorex_kwargs: A string of additional parameters for the `SKY` recipe from the MUSE Pipeline User Manual that are not included in the :obj:`json` configuration file otherwise can be provided
   :proptype esorex_kwargs: string
   :options esorex_kwargs: default='none'

   :property create_sof: Must be set :obj:`True` if one wants to create a new `.sof` file. In case the user wants to change or create the `.sof` file manually,
                        than it must be set to :obj:`False` since it will be overwritten otherwise
   :proptype create_sof: bool
   :options create_sof: false, true, default=true

.. json:object:: sci_post

   :property execute: Must be set :obj:`True` if this step should be executed
   :proptype execute: bool
   :options execute: false, true, default=true

   :property autocalib: This may execute the flux autocalibration between the different IFUs for empty fields. If ``user`` a AUTOCAL_FACTORS.fits table has to be provided by the user.
   :proptype autocalib: string
   :options autocalib: 'none', 'deepfield', 'user', default=none

   :property subtract_sky: The user can decide to not run any sky subtraction. All telluric lines will be included in the final datacube. this may be useful if one
                            wants to do their own wavelength calibration.
   :proptype subtract_sky: bool
   :options subtract_sky: false, true, default=true

   :property raman: If laser guide stars are used raman scattering in the atmosphere may be visible in the final data cubes. If set to :obj:`True` the Raman lines
                     are removed. For more details we refer to the `MUSE data reduction pipeline`_ handbook.
   :proptype raman: bool
   :options raman: false, true, default=false

   :property esorex_kwargs: A string of additional parameters for the `SCIPOST` recipe from the MUSE Pipeline User Manual that are not included in the :obj:`json` configuration file otherwise can be provided
   :proptype esorex_kwargs: string
   :options esorex_kwargs: default='none'

   :property create_sof: Must be set :obj:`True` if one wants to create a new `.sof` file. In case the user wants to change or create the `.sof` file manually it must be set to :obj:`False` since it will be overwritten otherwise
   :proptype create_sof: bool
   :options create_sof: false, true, default=true

.. json:object:: dither_collect

   :property execute: Must be set :obj:`True` if this step should be executed
   :proptype execute: bool
   :options execute: false, true, default=true

   :property user_list: A user list with the identifiers of the dither positions, in case only a subset shall be collected and copied to the final folder
                        to be used for the combination. If left empty all dither positions are used to create the final data cube.
   :proptype user_list: list
   :options user_list: default=[]

.. json:object:: exp_align

   :property execute: Must be set :obj:`True` if this step should be executed
   :proptype execute: bool
   :options execute: false, true, default=true

   :property esorex_kwargs: A string of additional parameters for the `EXP_ALIGN` recipe from the MUSE Pipeline User Manual that are not included in the :obj:`json` configuration file otherwise can be provided
   :proptype esorex_kwargs: string
   :options esorex_kwargs: default='none'

   :property create_sof: Must be set :obj:`True` if one wants to create a new `.sof` file. In case the user wants to change or create the `.sof` file manually,
                        than it must be set to :obj:`False` since it will be overwritten otherwise
   :proptype create_sof: bool
   :options create_sof: false, true, default=true

.. json:object:: exp_combine

   :property execute: Must be set :obj:`True` if this step should be executed
   :proptype execute: bool
   :options execute: false, true, default=true

   :property weight: The method how the fluxes in each dither position sre weighted. For more details we refer to
                     the `MUSE data reduction pipeline`_ handbook.
   :proptype weight: string
   :options weight: 'exptime', 'fwhm', 'none', 'header', default=exptime

   :property esorex_kwargs: A string of additional parameters for the `EXP_COMBINE` recipe from the MUSE Pipeline User Manual that are not included in the :obj:`json` configuration file otherwise can be provided
   :proptype esorex_kwargs: string
   :options esorex_kwargs: default='none'

   :property create_sof: Must be set :obj:`True` if one wants to create a new `.sof` file. In case the user wants to change or create the `.sof` file manually,
                        than it must be set to :obj:`False` since it will be overwritten otherwise
   :proptype create_sof: bool
   :options create_sof: false, true, default=true



Utility methods
---------------

The following methods are routines to support the :class:`MUSEreduce.musereduce`, needed to execute the various data reduction steps as it is suggested in the `MUSE data reduction pipeline`_. They are documented here for completness.

.. automodule:: MUSEreduce
   :members:
   :private-members:
   :exclude-members: musereduce
   :member-order: bysource

.. _Zeidler et al. 2019: https://ui.adsabs.harvard.edu/abs/2019AJ....158..201Z/abstract

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
   supports now *pipeline v2.4.2* and the *NFM-AO* added: pipeline_path choosing if darks may be used only reduces *STD* once per OB general use of external *SKY* fields collecting the files for ``muse_exp_combine`` in an independent step ``muse_exp_align`` is an independent step now
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
   :class:`MUSEreduce.musereduce` can now handle if the exposures for one pointing are distributed via multiple OBs with multiple exposures in each OB.
.. versionadded:: 1.0
   The release version as originally published in `Zeidler et al. 2019`_.
.. versionadded:: 1.0.2
   It is possible to omit the data reduction of the standard star in ``muse_scibasic`` by setting the *execute_std* to :obj:`False`.
   More skylines can be added in the ``muse_scibasic`` module but adding additional wavelengths in Angstrom to *skylines*.
.. versionadded:: 1.1.0
   :class:`MUSEreduce.musereduce` supports now *pipeline v2.8.1*. A legacy version for *pipeline v2.4.2* was created in a separate release. The keywords *autocalib* and in ``sci_oost`` was added.
   The naming convention in ``dither_collect`` was changed so that it matches the convention  **RADEC_EXPTIME_ROTANGLE_COUNTER**.
.. versionadded:: 1.1.1
   The key words *srcmin* and *srcmax* were added to *config.json* to be used in ``exp_align``.
.. versionadded:: 1.2.0
    :class:`MUSEreduce.musereduce` supports now *pipeline v2.8.7*. Additonally, kwargs can be set in the `json` configuration file for the different routines to support any other ``esorex`` command that is not specifically listed.
.. versionadded:: 1.2.1
    :class:`MUSEreduce.musereduce` changed that the sof file gets included individually when calling the sorex command. thats was needed that additional kwargs can be provided.