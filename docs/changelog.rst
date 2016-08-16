Changelog
=========

dev
---
 * Each bin formed by :meth:`~ifuanal.IFUCube.emission_line_bin` is now
   restricted to a contiguous group of pixels satisfying the algorithm
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
   - The first HDU is the primary extension of the cube. Only
     the header is read from this extension - any data, if it even exists, is
     not accessed.
   - The second HDU is the science data (NaNs for bad data). A WCS in the
     header of this extension is required for some functionality.
   - The third HDU is the standard deviation of the science data.
 * :ref:`Saving and loading instances <saving-loading>` has been updated to
   separate the cube data from the pickle file. The FITS file described by
   ``cube_hdu`` is saved separately from the pickle file with the extension
   `.fits`. This is reread when loading an instance.

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
