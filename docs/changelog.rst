Changelog
=========

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
