Changelog
=========

vdev
----
 * :attr:`results` dict has been reimplemented with a new structure and now
   supercedes :attr:`bin_nums`.

   - Each bin's spectrum and nucleus distances are now calculated upon bin
     creation.
   - Continuum and emission line fitting are contained in separate entries of
     each bin's results entry named ``"continuum"`` and ``"emission"``.
 * :meth:`parse_results` and :meth:`parse_emission` have been renamed to
   :meth:`~ifuanal.IFUCube._parse_continuum` and
   :meth:`~ifuanal.IFUCube._parse_emission`. They are automatically called
   after fitting so do not need to be explicitly called usually.
 * The Balmer decrement is calculated when parsing the emission line model and
   an estimate of E(B-V)_gas is made to deredden the fluxes of the emission
   lines.
 * Added :meth:`~ifuanal.IFUCube.plot_line_map` to visualise the results for a
   given emission line.
 * Added :meth:`~ifuanal.IFUCube.plot_bpt` to produce a BPT diagram.
 * Exposed the filepath of the json emission lines file via the ``el_json``
   argument to :class:`~ifuanal.IFUCube` and :class:`~ifuanal.MUSECube`.
 * Emission line fitting now fits the the line widths of each entry in
   `data/emission_lines.json` separately (doublet/triplet lines will be
   tied to the same widths). Previously only two widths were fit, one for
   forbidden lines and one for Balmer.
 * Custom bins are now highlighted on the radial and cumulative metallicity
   plots.
 * Statistical uncertainties on metallicity measurements are now calculated.
 * Emission line results for ``fwhm`` and ``offset`` (and their uncertainties)
   are now given in km/s.
 * Emission line results are now length 2 tuples of (value, uncert) and the
   defunct entries (e.g. ``flux_uncert``, ``fwhm_uncert``) are removed.
 * Initial bin seeds are now plotted on the bin map for emission line binning.
 * Fixed the behaviour of emission line bin to now reject bin seeds within
   ``border`` of a ``nan`` value rather than only within ``r`` of the edge of
   the array.
 * The header card ``BUNIT`` is now queried in the `DATA` extension to get the
   flux units in order to label plots.
 * The spectra returned by :meth:`_get_weighted_spectrum` are now multiplied
   by the number of spaxels in the bin to conserve flux.
 * Removed ``use_tmp_dir`` and ``append`` arguments from :meth:`run_starlight`.
 * Fixed bug in calculation of bin distances from nucleus.
 * Removed unused components in bundled bc03 base.
 * Added ``cumweight`` option to :meth:`plot_metallicity`. This will weight
   the cumulative metallicity plot by the H\ :math:`\alpha` flux of each bin,
   as a proxy for the SFR.
 * Added ``smooth`` parameter to :meth:`emission_line_bin` to smooth the line
   map before peak detection.
 * BPT plot is only made for bins with S/N > 3 for all relevant line fluxes.
 * Fix bug in bin spectra uncertainty calculation.
 * Added ``resid_sig`` argument to :meth:`run_emission_lines` to further remove
   broad continuum residuals.
 * Added N2 indicator from Marino+13 (``M13_N2`` in results dict) and renamed
   the existing ``M13`` (O3N2) indicator to ``M13_O3N2``.
 * Updated continuum level determination for EW measurements of emission lines
 * Added ``cont`` entry to emission line results, giving the continuum level
   used for the EW measurement.
 * Added new emission lines to STARLIGHT mask file.
 * Allow ``el_json`` argument to :class:`IFUCube` to now be given as a
   dictionary or json filename.
 * Apply Balmer decrement correction to the uncertainties of flux and
   equivalent width.

v0.3.0
------
 * Overhauled behind-the-scenes in emission line fitting to make more general
   for custom line additions:

   - The :class:`astropy.modeling.CompoundModel` is no longer stored due to
     issues with pickling, but is instead created as and when necessary via
     :func:`~ifuanal._get_emline_model`.
   - The emission line fitting results are all stored in the :attr:`results`
     dictionary entry designated by line name following
     `data/emission_lines.json`.
   - `data/emission_lines.json` can now be added to with other lines, although
     the existing entries should not be altered.
 * :mod:`pathos` is no longer a requirement.
 * Each bin formed by :meth:`~ifuanal.IFUCube.emission_line_bin` is now
   restricted to be a contiguous group of pixels satisfying the algorithm
   conditions. A minimum size of the bins in pixels can be set via the new
   ``min_npix`` argument.
 * Updated :meth:`~ifuanal.IFUCube.get_weighted_spectrum` to return the
   weighted mean of the uncertainty in the spaxels rather than the absolute
   differences from the mean.
 * Renamed :meth:`plot_kinematic` to :meth:`~ifuanal.IFUCube.plot_kinematics`
 * Added :meth:`~ifuanal.IFUCube.get_loc_bin`
 * Updated metallicity plotting to include more plots.
 * Fixed bug in O3N2 metallicity calculation.
 * :class:`~ifuanal.IFUCube` now takes ``cube_hdu`` as an argument, a
   :class:`astropy.io.fits.HDUList` (see `<here
   http://docs.astropy.org/en/stable/io/fits/api/hdulists.html>`_), instead of
   separate science data and stddev cube headers. ``cube_hdu`` should be
   length-3 with the following entries:

   1. Primary extension of the cube. Only the header is read from this
      extension - any data, if it even exists, is not accessed.
   2. Science data (NaNs for bad data). A WCS in the
      header of this extension is required for some functionality.
   3. Standard deviation of the science data.
 * :ref:`Saving and loading instances <saving-loading>` has been updated to
   separate the cube data from the pickle file. The FITS file described by
   ``cube_hdu`` is saved separately from the pickle file with the extension
   `.fits`. This is reread when loading an instance.
 * Merged :meth:`get_emission_line_metallicities` and
   :meth:`get_emission_line_fluxes` into
   :meth:`~ifuanal.IFUCube.parse_emission`.
 * Renamed :meth:`plot_emission_lines` into
   :meth:`~ifuanal.IFUCube.plot_emission`

v0.2.0
------
 * :ref:`HII region binning <hii-binning>` algorithm added.
 * Removed :attr:`vor_sn` as an attribute of :class:`IFUCube`. Replaced with
   ``target_sn`` argument to :meth:`~ifuanal.IFUCube.voronoi_bin`.
 * :meth:`plot_starlight_results()` renamed to
   :meth:`~ifuanal.IFUCube.plot_continuum`.
 * Reorganised :attr:`bin_num`. Now a dictionary of `mean` and `spax` entries
   for each bin. See :ref:`binning`.

v0.1.0
------
 * First release
