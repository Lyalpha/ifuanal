Tutorial
========

ifuanal requires a reduced (i.e. sky-subtracted, flux- and
wavelength-calibrated) datacube as ingestion. Wavelength scale should be in
angstroms.

This tutorial will use an example workflow that performs :ref:`stellar
continuum <cont-fitting>` and :ref:`emission line <emission-fitting>` fitting and
subsequent :ref:`analysis <analysis>` on a MUSE data cube.

A description of the processes in each step as well as some of the pertinent
arguments is given below. For a full description of optional arguments and
their format, see the :ref:`api documentation of the methods <ifuanal-api>` (or type
``object?`` from within an IPython session to find out about ``object``). All
pixel-coordinates should be passed to methods as zero-indexed in the order
``x``, ``y`` (in the FITS standard). Any files produced are defined in this
tutorial by their `suffix.extension`. The prefix will be that of the cube being analysed
(minus the filename extension).

.. NOTE::

   ifuanal has only been tested on MUSE data so far and the tutorial will follow
   the process for a MUSE data cube.

   Although in principle other instruments' data can be analysed in a similar
   fashion, a custom child class of :class:`~ifuanal.IFUCube` similar to
   :class:`~ifuanal.MUSECube` would need to be written to properly format the
   data to how :class:`~ifuanal.IFUCube` expects it.

.. NOTE::

   In this tutorial, the markup of different objects are as follows:
   
   * :class:`class_name`
   * :attr:`class_attribute`
   * :meth:`class_method`
   * ``variable_name`` (or ``argument_name``)


Create a new instance
---------------------

For this example MUSE science verification data of the target **NGC2906** will be analysed.

::

  >>> import ifuanal 
  >>> cube = ifuanal.MUSECube("NGC2906.fits", 0.008138)

This will initialise the :class:`~ifuanal.MUSECube` class, which does some
small manipulation to the MUSE FITS file input before ingestion to
:class:`~ifuanal.IFUCube`, namely:

* Split the MUSE FITS file into its `DATA` and `STAT` extensions.
* Add a header card `IFU_EBV` specifying the reddening. The argument ``ebv`` can
  be passed to :class:`~ifuanal.MUSECube` to explicitly set this, otherwise its default value
  of "IRSA" will contact the Infrared Science Archive to automatically determine
  it based on the coordinates of the WCS reference pixel of the cube (this
  requires the optional dependancy :mod:`astroquery` to be installed).
* Add a header card `IFU_Z` specifying the redshift. In the example case this is
  `0.008138`
* The MUSE data `STAT` extension gives the variance of the science
  data. ``IFUCube`` wants the standard deviation and so we square root this
  extension.

``IFUCube`` is then initialised which will set up the wavelength scale, check
the STARLIGHT directory (:attr:`sl_dir`) exists, and load the emission line data
from `data/emission_lines.json`

.. WARNING::

   The lines in `data/emission_lines.json` should not be altered in the current
   version. There are some hard coded sections that require these and only these
   lines. A generalising of this process will allow custom line lists to be
   given.

.. _deredden-deredshift:

Deredden and deredshift
-----------------------

::

  >>> cube.deredden()
  dereddening with E(B-V) = 0.040mag and RV = 3.1
  >>> cube.deredshift()
  deredshifting from z = 0.008138

These are pretty self-explanatory. One thing to note is that the `E(B-V)` and
`z` values are taken from header cards ``IFU_EBV`` and ``IFU_Z``,
respectively. Dereddening is done using a Cardelli, Clayton and Mathis (1989)
polynomial.

Once either method has been called the appropriate header values is
set to `0` and subsequent calls will not do anything to the cube, e.g.::

  >>> cube.deredden()
  ebv = 0, skipping deredden()

The wavelength array attribute :attr:`lamb` is updated with the deredshifting:
::

  >>> print("{:.2f}, {:.2f}".format(cube.lamb[0], cube.lamb[-1]))
  4711.66, 9274.52

Mask foregound/background sources
---------------------------------

We can remove spaxels from the data cube (by setting their values to ``np.nan``)
to ensure they are not considered in subsequent analysis. For NGC2906 there is a
foreground star in our cube, which we want to mask: ::

  >>> cube.mask_regions([(109, 192),], 12)
  masking regions

``109, 192`` are the approximate pixel coordinates of the star
and ``12`` is the radius of the mask in pixels. Note the coordinates of the
regions should be given as a list of length-2 lists/tuples. The radius argument
can be a list also, in order to specify a different radius for each region to
mask, or, if ``len(regions) > len(radii)`` it will loop over the radii. e.g. for
multiple regions: ::

  >>> # cube.mask_regions([(10, 20), (30, 40), (50, 60)], [8, 9, 10])

will use radii of ``8``, ``9`` and ``10`` for the three regions, whereas: ::

  >>> # cube.mask_regions([(10, 20), (30, 40), (50, 60)], 10)

will use a radius of ``10`` for all regions.


Find the galaxy centre
----------------------

We need to provide an initial guess to find centre of the galaxy, usually by
simply eyeballing the cube. This can be given in pixel coordinates or RA and
DEC if the argument ``usewcs = True``. The centre is found by fitting a 2D
gaussian to a region around this initial guess.

To resolve poor fits look at the docs for :meth:`~ifuanal.IFUCube.set_nucleus`,
since there are other arguments to play with, as well as the option to specify
a location outside the FOV. ::

  >>> cube.set_nucleus([162, 167])
  set nucleus as (160.592, 166.442)

By default this will also produce a plot `_nucleus.pdf` showing the data, model
and residual for checking (``plot=False`` to skip this).

.. TODO::

   The use of this in the analsis is currently quite limited. Further updates
   will use this to calculate e.g. deprojected distances of bins and provide
   maps in terms of offset from the centre.

Binning the spaxels
-------------------

We do not want to consider sky spaxels in our analysis and, additionally, we do
not want to perform fitting to low signal-to-noise ratio (SNR) spaxels. To
circumvent this we employ spaxel binning.

The spaxels are to be binned into distinct regions in order to increase the S/N
of the composite region spectra for fitting. :ref:`vor-binning` and
:ref:`hii-binning` are the two currently implemented methods, with the ability
to also :ref:`add custom bins <custom-bins>`.

These binning routines will populate an attribute of :class:`~ifuanal.IFUCube`
named :attr:`bin_nums`. Once a binning routine has been performed, this will
appear as a ``dict`` of the form: ::

  >>> cube.bin_nums
  {0: {'mean': (x_mean, y_mean),
  'spax': (x_spax, y_spax)},
   1 ...
  }

The entry ``mean`` gives the centre of the bin (for Vornoi binning this is the
centre of mass, for the HII region binning, this is the peak). ``spax`` gives
the spatial indices of the datacube belonging to that bin. e.g. ::

  >>> cube.bin_nums[100]["spax"] # get the x and y coordinates of bin 100

.. Note::

   To repeat or redo binning once :attr:`bin_nums` is populated, pass the
   argument ``clobber= True`` in the binning method's call.

.. _vor-binning:

Voronoi binning
^^^^^^^^^^^^^^^

`Voronoi tessellation <https://en.wikipedia.org/wiki/Voronoi_diagram>`_ is
performed using the `Voronoi binning algorithm
<http://www-astro.physics.ox.ac.uk/~mxc/software/>`_ to produce bins from
spaxels with individual S/N > 3. The individual spectra in each bin are
combined to increase the SNR to some target value.

The SNR of the spectra are calculated in a specific wavelength window (default
is 5590 to 5680) and emission line signal-to-noise ratios can be estimated by
subtracting off a continuum SNR (see docs for
:meth:`~ifuanal.IFUCube.voronoi_bin`) ::

  >>> cube.voronoi_bin()
  voronoi_binning with S/N target of 20
  ... [voronoi output] ...


By default a plot of the bins and their S/N will be produced and saved as
`_bins_vor.pdf`.

.. _hii-binning:

HII region binning
^^^^^^^^^^^^^^^^^^

This binning algorithm uses the method of `HII explorer
<http://www.caha.es/sanchez/HII_explorer/>`_, with a python
implementation, to grow bins around peaks in the emission line flux. ::

  >>> cube.emission_line_bin(min_peak_flux=1100, min_frac_flux=0.1,
  ... max_radius=5, min_flux=600)

A description of these required arguments is available in the documentation for
:meth:`~ifuanal.IFUCube.emission_line_bin`. These will have to be tailored to 
each data cube. Although usually the binning will be done for the
H\ :math:`\alpha` line, any line or wavelength can be chosen via the
``line_lamb`` argument.

Briefly, the method is:

1. :func:`~ifuanal.get_line_map` is called. This returns an emission line map
by simulating a narrowband filter observation of the datacube and subtracting a
continuum determined by two neighbouring filters.

2. All peaks equal to or above ``min_peak_flux`` in this map are found via
:func:`scipy.ndimage.maximum_filter`. These peaks are allowed to be close as
the subsequent growth of the bins will merge nearby peaks.

3. Starting with the brightest, these peaks are the seeds for new bins. All
   nearby pixels that satisfying the following are included in the bin:

   * within ``max_radius`` of peak.
   * flux is above ``min_flux`` and ``min_frac_flux`` :math:`\times` peak
     flux.
   * fas not already been allocated a bin.

The resulting bins are then saved in :attr:``bin_nums``. By default a plot of
the emission line map creation and the bins will be produced and saved as
`_bins_el.pdf`.

.. _custom-bins:

Adding custom bins
^^^^^^^^^^^^^^^^^^

Custom bins can be added by defining a centre and radius. These bins will have
negative bin numbers beginning at ``-1`` in results.

As an example we make an SDSS-like 3 arcsec fibre on the galaxy nucleus: ::

  >>> cube.add_custom_bin([160.592, 166.442], 3/0.2)
  "added custom bin -1 to the list"

where 0.2 is the pixel scale of MUSE in arcsecs. Once all fitting has been performed, the results for this bin (assuming it was the first or only custom bin to be added) can be accessed via the bin number -1 in the :ref:`results-dict`

.. TODO::

   Currently this is limited only to circular bins but arbitrary bins (by just
   specifying a list of ``x`` and ``y`` pixel coordinates) should be added.

.. WARNING::

   Where spaxels are included in multiple bins, the plots will not represent
   these correctly (or consistently?).

.. _cont-fitting:

Stellar continuum fitting
-------------------------

Stellar continuum fitting is performed via `STARLIGHT
<http://astro.ufsc.br/starlight/>`_ (see :ref:`starlight-install`). 

**The tl;dr version:** ::

  >>> cube.run_starlight()
  running starlight fitting
  ... [starlight fitting output] ...
  >>> cube.parse_results()
  ... [parsing output] ...

**Extended version:**

Recommended reading for more information on the setup of STARLIGHT and in
particular the format of the config/mask/grid files is the extensive manual for
version 4.

By default all bins will be fitted, or a list of bin numbers can be passed
explicitly as the :attr:`bin_num` argument. The default set of bases are 45
Bruzual & Charlot (2003) models, this can be changed through the use of the
``base_name`` argument and the inclusion of the appropriate files in
:attr:`sl_dir` (see below). A temporary directory is also created
`/tmp/starlight_[random]` to store all the output.

:meth:`~ifuanal.IFUCube.run_starlight` searches :attr:`sl_dir` (default is
`starlight/` subdir of ifuanal\'s directory) for the following files:

* `starlight.config` - the main configuration file for the STARLIGHT
  run. In particular it contains limits on fittable values and specifies the
  wavelength window for normalisation of the spectra. The default config file
  with ifuanal is set up for a balance of robust fitting and speed.
* `starlight.mask` - a list of wavelength windows (around emission lines) to
  mask in the fitting of the continuum.
* a directory named ``base_name`` and a file named '``base_name``\ `.base`' -
  the choice of base models to use as well as the directory containing the bases
  (both must exist with these naming formats for ``base_name`` to be valid). We
  resample the bases to the same wavelength step as our deredshifted data cube
  (to avoid manipulating our data and introducing correlated uncertainties).

The process for a single bin is as follows:

1. Find all spaxels in the bin and
   a. if we have 1 spaxel, return that spectrum from the data cube.
   b. if we have >1 spaxels, return the weighted mean of these spectra from the
   data cube.
2. Write the spectrum to `/tmp/starlight_[random]/spec_[random]`.
3. Write a `grid` file used by STARLIGHT to
   `/tmp/starlight_[random]/grid_[random]`.
4. Call the STARLIGHT executable for this bin and return the file name of the
   output.

Once all bins are fit the file names of the results are mapped to the bin
numbers in a text file `_sloutputs.txt` as well as in a dictionary,
:attr:`sl_output`, in the class.

A call to :meth:`~ifuanal.IFUCube.parse_results` will then read these STARLIGHT output files
and store the pertinant results into a large dictionary :attr:`results` (see
:ref:`results-dict`). Any bins without output or where the output does not
follow the standard STARLIGHT output style will not be saved and their values
will also be removed from :attr:`bin_nums` (i.e. they are not considered in the
emission line fitting).

.. _emission-fitting:

Emission line fitting
---------------------

Emission line fitting is done with a set of single gaussians, one for each of
the lines given in `data/emission_lines.json`.

**The tl;dr version:** ::

  >>> cube.run_emission_lines()
  ... [emission lines output] ...
  >>> cube.get_emission_line_fluxes()
  >>> cube.get_emission_line_metallicities()

**Extended version:**

The emission line model is formed from the addition of gaussians via
`astropy\'s compound models
<http://docs.astropy.org/en/stable/modeling/compound-models.html>`_ and is fit
using a `Levenberg-Marquardt LSQ fitter
<http://docs.astropy.org/en/stable/api/astropy.modeling.fitting.LevMarLSQFitter.html#astropy.modeling.fitting.LevMarLSQFitter>`_.

As with the :ref:`continuum fitting <cont-fitting>`, by default all bins (that
have a valid STARLIGHT output) are fit, or a list of specific bins to be fit
can be passed as ``bin_num``.

Especially with lower SNR features, the fitter is suceptible to finding local
minima in the LSQ sense and is very sensitive to the inital guess for the
amplitude, mean and standard deviation of the gaussians. To circumvent this a
somewhat brute force method is overlaid on the minimisation of the fitter, as
well as applying some conditions to the fitted parameters:

* A grid of initial guesses with every combination of the initial guess lists is formed:
  * The arguments ``vd_init``, ``v0_init`` and ``amp_init`` are the
  initial guesses for the standard deviation (in km/s), mean offset (in km/s)
  and amplitude (in units of ``fobs_norm`` -- see STARLIGHT). See the docs
  for :meth:`~ifuanal.IFUCube.run_emission_lines` for more information.
* The standard deviation of the emission lines are restricted to between 5 and
  120 km/s by default, this can be altered with the argument
  ``stddev_bounds``.
* The offset of the lines is limited to between -500 and +500 km/s (from the
  overall :ref:`deredshifted cube <deredden-deredshift>`) by
  default, this can be altered with the argument ``offset_bounds``.
* The mean and standard deviation of the balmer lines are tied to be the
  same. The forbidden lines are also tied to each other but they can differ
  from the balmer values.
* If any negative amplitude is found, it is set to zero (since we are dealing
  only with emission lines currently).

Each of the initial guess combinations in the grid is fitted with the fitter and the :math:`\chi^2`/dof value of the fit stored; the minimum :math:`\chi^2`/dof is taken as the best fit.

.. _saving-loading:

Saving and loading instances
----------------------------

It is possible to save your current instance to preserve results and load these
results later or elsewhere via pickling (performed with `dill
<`https://github.com/uqfoundation/dill>`_). ::

  >>> cube.save_pkl()
  writing to temporary pickle file /cwd/ifuanal_[random].pkl
  moving to NGC2906.pkl

The instance ``cube`` is now stored in `NGC2906.pkl`, including all results of
fitting, the data and stddev cube etc. This file can then be loaded later to
resturn to the same state:::

  >>> cube2 = ifuanal.IFUCube.load_pkl("NGC2906.pkl")
  loaded pkl file NGC2906.pkl

And ``cube2`` will have all the attributes of the ``cube`` class: ::

  >>> print cube2.nucleus
  (160.592, 166.442)

.. WARNING::

   This appears to breakdown with large numbers of bins (not sure exactly but
   >6000). This is probably a filesize issue that may be fixed by saving the
   altered datacube separately from the pickle file and reingesting it upon load
   since.


.. _analysis:

Analysis
--------

After the fitting has been done for the continuum and emission lines, then we
can do all this fancy stuff...

.. _results-dict:

:attr:`results` dictionary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
As an example, to see the results for a bin of number
``bn``, type:::

  >>> cube.results["bin"][bn]
  ... [bin results output] ...

The ``results`` dictionary contains

.. TODO::
   Write this section.


Plotting
^^^^^^^^

Once all fitting has been done, maps of the bins with various quantities can be
output...

.. TODO::
   Write this section.

.. Warning:: 
   Plotting will fail if use more than 6 metallicities are used in the
   STARLIGHT bases, or if the number of ages for each metallicity are different.
