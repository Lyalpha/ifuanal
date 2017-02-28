Tutorial
========

ifuanal requires a reduced (i.e. sky-subtracted, flux- and
wavelength-calibrated) datacube as ingestion.

This tutorial will use an example workflow that performs :ref:`stellar
continuum <cont-fitting>` and :ref:`emission line <emission-fitting>` fitting
and subsequent :ref:`analysis <analysis>` on a MUSE data cube.

A description of the processes in each step as well as some of the pertinent
arguments is given below. For a full description of optional arguments and
their format, see the :ref:`api documentation of the methods <ifuanal-api>` (or
type ``object?`` from within an IPython session to find out about
``object``).

* All pixel-coordinates should be passed to methods as zero-indexed in the order
  ``x``, ``y`` (in the FITS standard).
* Any files produced are defined in this tutorial by their
  `suffix.extension`. The prefix will be that of the cube being analysed (minus
  the filename extension).
* By default, ifuanal will spawn ``min(Ncpu-1, Ncpu*0.9)`` processes for the
  multiprocessed parts of the analyses where ``Ncpu`` is the number of cpus on
  the system (this is saved under the attribute :attr:`n_cpu` and can be changed
  manually).
* Emission lines to fit are defined in the file `data/emissionlines.json`.
  Additional lines can be fit for by copying this file and manually adding new
  entries, and giving the filepath of this new file via the ``el_json``
  argument (alternatively a pythin ``dict`` can be passed of the same format)
  **If the existing entries in the default line list are removed or altered
  some plotting may fail that relies on these lines**. No sanity checking of
  the lines is done, so make sure they are sensible for the data being analysed
  (i.e. within the wavelength window of the cube).

.. NOTE::

   ifuanal has only been tested on MUSE data so far and the tutorial will follow
   the process for a MUSE data cube.

   Although in principle other instruments' data can be analysed in a similar
   fashion, a custom child class of :class:`~ifuanal.IFUCube` similar to
   :class:`~ifuanal.MUSECube` would need to be written to properly format the
   data to how :class:`~ifuanal.IFUCube` expects it.

.. NOTE::

   In this tutorial, the markup of different objects are as follows:

   * :class:`Class_Name`
   * :attr:`class_attribute`
   * :meth:`class_method`
   * ``variable_name`` (or ``argument_name``)


Create a new instance
---------------------

For this example MUSE science verification data of the target **NGC2906** will
be analysed.

::

  >>> import ifuanal
  >>> cube = ifuanal.MUSECube("NGC2906.fits", 0.007138)

This will initialise the :class:`~ifuanal.MUSECube` class, which does some
small manipulation to the MUSE FITS file input before ingestion to
:class:`~ifuanal.IFUCube`, namely:

* Open the MUSE FITS file into a :class:`astropy.io.fits.HDUList` of the
  `PRIMARY`, `DATA` and `STAT` extensions.
* Add a `PRIMARY` header card `IFU_EBV` specifying the reddening. The argument
  ``ebv`` can be passed to :class:`~ifuanal.MUSECube` to explicitly set this,
  otherwise its default value of "IRSA" will contact the Infrared Science
  Archive to automatically determine it based on the coordinates of the WCS
  reference pixel of the cube (this requires the optional dependancy
  :mod:`astroquery` to be installed).
* Add a PRIMARY header card `IFU_Z` specifying the redshift. In the example
  case this is `0.007138`
* The MUSE data `STAT` extension gives the variance of the science
  data. ``IFUCube`` wants the standard deviation and so we square root this
  extension.

``IFUCube`` is then initialised - this will set up the wavelength scale, check
the STARLIGHT directory (:attr:`sl_dir`) exists, and load the emission line data
from ``el_json`` (default `data/emission_lines.json`).

.. NOTE::

   The input FITS file must contain the header cards ``CUNIT3`` and ``BUNIT``
   in the `DATA` extension, which are parsable by :mod:`astropy.units`\' `string
   parser
   <http://docs.astropy.org/en/stable/units/format.html#creating-units-from-strings>`_.

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

To correct bad fits, look at the docs for :meth:`~ifuanal.IFUCube.set_nucleus`,
since there are other arguments to play with, as well as the option to specify
a location outside the FOV. ::

  >>> cube.set_nucleus(162, 167)
  set nucleus as (160.592, 166.442)

By default this will also produce a plot `_nucleus.pdf` showing the data, model
and residual for checking (``plot=False`` to skip this).

.. TODO::

   The use of this in the analsis is currently quite limited. Further updates
   will use this to calculate e.g. deprojected distances of bins and provide
   maps in terms of offset from the centre.

.. _binning:

Binning the spaxels
-------------------

We do not want to consider sky spaxels in our analysis and, additionally, we do
not want to perform fitting to low signal-to-noise ratio (SNR) spaxels. To
circumvent this we employ spaxel binning.

The spaxels are to be binned into distinct regions in order to increase the S/N
of the composite region spectra for fitting. :ref:`hii-binning` and
:ref:`vor-binning` are the two currently implemented methods, with the ability
to also :ref:`add custom bins <custom-bins>`.

These binning routines will populate :ref:`results-dict` with each bin. The
information is stored as follows for bin number ``bn``: ::

  >>> cube.results["bin"][bn]
  {'mean': (x_mean, y_mean),  # the pixel coordinates of the centre of the bin
   'spax': (x_spax, y_spax)}, # the pixel coordinates of the spaxels in the bin
   'spec': 4xN array,         # cols: lambda, flux(lambda), sigma(lambda), flag
   'dist_min': float,         # minimum distance to nucleus
   'dist_max': float,         # maximum distance to nucleus
   'dist_mean': float,        # distance of 'mean' to nucleus
   'continuum': {},           # dict populated once continuum fitting is done
   'emission': {},            # dict populated once emission fitting is done
  }

For Vornoi binning, ``mean`` is the centre of mass, whereas for the HII region
binning, this is the seed peak.

In the case of a single spaxel bin, ``spec`` is just copied from the input data
and stddev cube. For a multi-spaxel bin, the weighted mean of the data and
uncertainties of all individual spaxels in the bin are used.

See :ref:`results-dict` for information on accessing and using this information.

.. Note::

   To repeat or redo binning, pass the argument ``clobber= True`` in the
   binning method's call. **This will also remove existing bin results
   including continuum and emission fitting.**

.. _hii-binning:

HII region binning
^^^^^^^^^^^^^^^^^^

This binning algorithm uses the method of `HII explorer
<http://www.caha.es/sanchez/HII_explorer/>`_, with a python
implementation, to grow bins around peaks in the emission line flux. ::

  >>> cube.emission_line_bin(min_peak_flux=1100, min_frac_flux=0.1,
  ... max_radius=5, min_flux=600)
  binning spaxels using HII explorer algorithm around emission line 6562.8
  processing bin seed [i]/[m]
  found [n] bins

A description of these required arguments is available in the documentation for
:meth:`~ifuanal.IFUCube.emission_line_bin`. These will have to be tailored to
each data cube. Although usually (and by default) the binning will be done for
the H\ :math:`\alpha` line, any line or wavelength can be chosen via the
``line_lamb`` argument.

Briefly, the method is:

1. :func:`~ifuanal.get_line_map` is called. This returns an emission line map
by simulating a narrowband filter observation of the datacube and subtracting a
continuum determined by two neighbouring filters.

2. The emission line map is filtered with a gaussian, whose width is given by
   the ``smooth`` argument. This helps to avoid picking noise peaks in the
   wings of bright emission regions, but this can be skipped by setting
   ``smooth`` to zero.

3. All peaks equal to or above ``min_peak_flux`` in the emission line map are
found via :func:`scipy.ndimage.maximum_filter`. These peaks are allowed to be
close since the subsequent growth of the bins will merge nearby peaks.

4. Starting with the brightest, these peaks are the seeds for new bins. All
   nearby pixels that satisfying the following are included in the bin:

   * within ``max_radius`` of peak.
   * flux is above ``min_flux`` and ``min_frac_flux`` :math:`\times` peak
     flux.
   * is not already been allocated a bin.

The resulting bins are then saved in ``cube.results["bin"]``. By default a
plot of the emission line map creation and the bins will be produced and saved
as `_bins_el.pdf`.

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

  >>> cube.voronoi_bin(target_sn=20)
  binning spaxels with Voronoi algorithm with S/N target of 20
  [voronoi output]
  processing bin [i]/[n]
  found [n] bins


The resulting bins are then saved in ``cube.results["bin"]``. By default a
plot of the bins and their S/N will be produced and saved as `_bins_vor.pdf`.

.. _custom-bins:

Adding custom bins
^^^^^^^^^^^^^^^^^^

Custom bins can be added by defining a centre and radius. These bins will have
negative bin numbers beginning at ``-1`` in results.

As an example we make an SDSS-like 2 arcsec fibre on the galaxy nucleus: ::

  >>> cube.add_custom_bin([160.592, 166.442], 2/0.2)
  "added custom bin -1 to the list"

where 0.2 is the pixel scale of MUSE in arcsecs. Once all fitting has been
performed, the results for this bin (assuming it was the first custom bin to be
added) can be accessed via the bin number -1 in the :ref:`results-dict`

.. TODO::

   Currently this is limited only to circular bins but arbitrary bins (by just
   specifying a list of ``x`` and ``y`` pixel coordinates) should be added.

.. WARNING::

   Where spaxels are included in multiple bins, the 2D map plots will not
   represent these correctly (or consistently?).

.. _cont-fitting:

Stellar continuum fitting
-------------------------

Stellar continuum fitting is performed via `STARLIGHT
<http://astro.ufsc.br/starlight/>`_ (see :ref:`starlight-install`).

**The tl;dr version:** ::

  >>> cube.run_starlight()
  running starlight fitting
  fitting [n] bins...
  STARLIGHT tmp directory for this run is /tmp/starlight_[random]/
  resampling base files [i]/[m]
  fitting bin number [i]
  parsing results
  [failed to parse /tmp/starlight_[random]/spec_[random]_out for bin [j]]
  parsing starlight output [i]/[n]

**Extended version:**

Recommended reading for more information on the setup of STARLIGHT and in
particular the format of the config/mask/grid files is the extensive manual for
version 4 `here <http://www.starlight.ufsc.br/papers/Manual_StCv04.pdf>`_.

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

1. Access the spectrum of the bin via :ref:`results-dict`.
2. Write this spectrum to `/tmp/starlight_[random]/spec_[random]`.
3. Write a `grid` file used by STARLIGHT to
   `/tmp/starlight_[random]/grid_[random]`.
4. Call the STARLIGHT executable for this bin and return the file name of the
   output (the spectrum file with a `_out` suffix).

Once all bins are fit, a call to :meth:`~ifuanal.IFUCube._parse_continuum` then
reads these STARLIGHT output files and parses the information into the
`"continuum"` entry in :attr:``results`` for each bin (see
:ref:`results-dict`). The dictionary entry `"continuum"` is populated with the
results of the STARLIGH fitting, please consulte the STARLIGHT documentation
(section 6 of the version 4 manual) for more information on these. In
particular, `"bases"` is the population mixture of the bases used to create the
best fitting continuum and `"sl_spec"` is the synthetic spectrum.

Any bins without output or where the output does not follow the standard
STARLIGHT output style will be shown in the terminal (`failed to
parse...`). This is usually due to normalisation errors in STARLIGHT where
there is ~0 flux in the continuum - the file printed to the terminal can be
inspected for further investigation. For a failed bin number of ``bn``, the
follow flag is set: ::

  >>> cube.results["bin"][bn]["continuum"]["bad"]
  1

This is ``0`` otherwise.



.. _emission-fitting:

Emission line fitting
---------------------

Emission line fitting is done with a set of single gaussians, one for each of
the lines given in ``el_json`` (default `data/emission_lines.json`).

**The tl;dr version:** ::

  >>> cube.run_emission_lines()
  fitting emission lines to [n] bins...
  fitting bin number [i]
  [no covariance matrix computed for bin [j], cannot compute fit uncertainties]
  emission line fitting complete
  parsing emission model [i]/[n]

**Extended version:**

The emission line model is formed from the addition of gaussians via
`astropy\'s compound models
<http://docs.astropy.org/en/stable/modeling/compound-models.html>`_ and is fit
using a `Levenberg-Marquardt LSQ fitter
<http://docs.astropy.org/en/stable/api/astropy.modeling.fitting.LevMarLSQFitter.html#astropy.modeling.fitting.LevMarLSQFitter>`_.

As with the :ref:`continuum fitting <cont-fitting>`, by default all bins (that
have a valid STARLIGHT output) are fit, or a list of specific bins to be fit
can be passed as ``bin_num``.

Especially with lower SNR features, the fitter is susceptible to finding local
minima in the LSQ sense and is sensitive to the inital guess for the
amplitude, mean and standard deviation of the gaussians. To circumvent this a
somewhat brute force method is overlaid on the fitter minimisation, as
well as applying some conditions to the fitted parameters:

* The residual spectrum is constructed by subtracting the continuum fit from
  the observed spectrum. This is then convolved with a gaussian of width
  ``resid_sig`` if required to further remove broad residuals (see below).
* The residual spectrum is masked for wavelengths more than ``offset_bounds`` +
  3 :math:`\times` ``stddev_bounds`` from an emission line rest wavelength.
  Wavelengths outside these windows are not fit for.
* A grid of initial guesses with every combination of the initial guess lists
  is formed. The arguments ``vd_init``, ``v0_init`` and ``amp_init`` are the
  initial guesses for the standard deviation (in km/s), mean offset (in km/s)
  and amplitude (in units of ``fobs_norm`` -- see STARLIGHT). See the docs
  for :meth:`~ifuanal.IFUCube.run_emission_lines` for more information.
* The standard deviation of the emission lines are restricted to between 5 and
  120 km/s by default, this can be altered with the argument
  ``stddev_bounds``.
* The offset of the lines is limited to between -500 and +500 km/s (from the
  overall :ref:`deredshifted cube <deredden-deredshift>`) by
  default, this can be altered with the argument ``offset_bounds``.
* The offset of the balmer lines are tied to be the same. The forbidden lines
  are also tied to each other but they can differ from the balmer values.
* The standard deviation width of the fits can differ between lines, but any
  doublets (or triplets) are forced to be fit with the same width.
* If any negative amplitude is found, it is set to zero (since we are dealing
  only with emission lines currently).

Each of the initial guess combinations in the grid is fitted with the fitter
and the :math:`\chi^2`/dof value of the fit stored; the minimum
:math:`\chi^2`/dof is taken as the best fit.

Parameters and their uncertainties are stored within the :ref:`results-dict`.
When a fitting is deemed to fail (``no covariance matrix computed for bin i
...``) this is either due to an inherently low SNR emission line spectrum or
the fitting encounted one of the bounding conditions of the fit
(``stddev_bounds`` or ``offset_bounds``). In the latter case the covariance of
the fitted parameters cannot be computed and an inspection of the fit via: ::

    >>> cube.plot_emission(i)
    plot saved to NGC2906_el_fit_i.png

will show the issue. For example, if the lines are well offset in velocity from
the galaxy, relaxing ``offset_bounds`` and providing ``v0_init`` with more
appropriate initial guesses should help the fit.

In the case of broad continuum residuals that are affecting the fitting,
these can be removed somewhat arbitrarily by using the argument ``resid_sig``.
This sets the sigma width in wavelength units of a gaussian kernel which the
residual emission line spectrum is convolved with. This convolved function
is then removed from the spectrum prior to fitting. It is important to not
fit the emission line of interest with this convolved function so ``resid_sig``
should be much larger than their widths.

.. _saving-loading:

Saving and loading instances
----------------------------

It is possible to save your current instance to preserve results and load these
results later or elsewhere via pickling (performed with `dill
<`https://github.com/uqfoundation/dill>`_). ::

  >>> cube.save_pkl()
  writing cube to temporary file /cwd/ifuanal_[random].pkl.fits
  moving to NGC2906.pkl.fits
  writing to temporary pickle file /cwd/ifuanal_[random].pkl
  moving to NGC2906.pkl

The instance ``cube`` is now stored in `NGC2906.pkl`, including all results of
fitting etc. Since problems can occur with very large pickle files, the cube
data is stored separately as a FITS file with the extension `.fits` added to
the pickle filename. This is a dereddened, deredshifted copy of the original
FITS file we loaded. A FITS file with the specific name `[pkl_filename].fits`
will be searched for when loading the instance and so a copy should be left
alongside the pickle file.

The instance can then be loaded later to return to the same state, by
specifiying the pickle file to load:::

  >>> cube2 = ifuanal.IFUCube.load_pkl("NGC2906.pkl")
  loaded pkl file NGC2906.pkl

And ``cube2`` will have all the attributes of the ``cube`` class, e.g.: ::

  >>> print cube2.nucleus
  (160.592, 166.442)

.. NOTE::

   The attribute :attr:`n_cpu` is updated upon loading an instance to be
   appropriate for the system being used.

.. _analysis:

Analysis
--------

After the fitting has been done for the continuum and emission lines, then we
can do all this fancy stuff...

.. _results-dict:

the :attr:`results` dictionary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
As an example, to see the results for a bin of number
``bn``, type: ::

  >>> cube.results["bin"][bn]

The ``results`` dictionary contains...

.. TODO::
   Write this section.


plotting
^^^^^^^^

Once all fitting has been done, maps of the bins and the results of the fitting
methods can be made.See the docs for each method for more info.

:meth:`~ifuanal.IFUCube.plot_continuum`
"""""""""""""""""""""""""""""""""""""""
plots the spectra of a bin and the fit to the continuum, as well as the
contribution of the various age and metallicity bases to the integrated fit.

:meth:`~ifuanal.IFUCube.plot_emission`
"""""""""""""""""""""""""""""""""""""""
plots the spectra of a bin and the fit to the emission spectrum.

:meth:`~ifuanal.IFUCube.plot_worst_fits`
""""""""""""""""""""""""""""""""""""""""
plots the ``N`` worst fits of each of the continuum and emission fits, as
determined by their :math:`\chi^2`/dof.

:meth:`~ifuanal.IFUCube.plot_yio`
"""""""""""""""""""""""""""""""""
plots the contribution of young, intermediate, and old stellar populations to
the continuum fits as a map.

:meth:`~ifuanal.IFUCube.plot_kinematics`
""""""""""""""""""""""""""""""""""""""""
plots the velocity offset and dispersion of the stellar populations in the
the continuum fits as a map.

:meth:`~ifuanal.IFUCube.plot_metallicity`
"""""""""""""""""""""""""""""""""""""""""
plots the metallicity for the chosen indicator as a map aslongside the
cumulative metallicity of the bins and a radial dependancy plot. If the
argument ``cumweight`` is ``True`` then the cumulative plot is weighted
by the SFR of each bin (i.e. the H\ :math:`\alpha` flux).
Custom bins are highlighted.

:meth:`~ifuanal.IFUCube.plot_line_map`
""""""""""""""""""""""""""""""""""""""

plots the EW, flux, velocity offset and FWHM of the chosen line (see
:meth:`~ifuanal.IFUCube.plot_line_map` docs).

:meth:`~ifuanal.IFUCube.plot_bpt`
"""""""""""""""""""""""""""""""""
plots the BPT diagram for each bind along with a 2D map of bin classifications.
Classification lines are taken from Kewley et al. (2013, ApJL, 774, 10) for the
AGN-HII division and Kewley et al. (2001, ApJ, 556, 121) for the maximal star
burst.

.. Warning::
   :meth:`~ifuanal.IFUCube.plot_continuum` will fail if use more than
   6 metallicities are used in the STARLIGHT bases, or if the number of ages
   for each metallicity are different.
