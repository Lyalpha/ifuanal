"""
IFUANAL

For the analysis of IFU data cubes.
"""

from __future__ import print_function

__version__ = "0.4.dev"
__author__ = "J. Lyman"

from itertools import repeat, cycle, product
import json
import math
import operator
import os
import re
import tempfile
import shutil
import subprocess
import sys
import warnings

from astropy.constants import c
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.modeling import models, fitting
import astropy.units as u
import astropy.wcs as wcs
import matplotlib.pyplot as plt
from matplotlib import gridspec, ticker, cm, colors, rc
import numpy as np
from scipy.interpolate import interp1d
from scipy import ndimage
from scipy.spatial.distance import cdist

from voronoi import voronoi
import dill as pickle
import multiprocessing as mp

# Stop numpy RuntimeWarnings (mainly to do with nans)
np.seterr("ignore")

# The file and directory path to this script
FILEPATH = os.path.realpath(__file__)
FILEDIR = os.path.dirname(FILEPATH)

# Colours for plots
OBSCOL = "k"  # observed spectrum
SYNTHCOL = "#D95F02"  # starlight fit
MASKCOL = "#7570B3"  # masked regions in fit
CLIPPEDCOL = "red"  # clipped by starlight
ERRCOL = "#666666"  # stddev
ZCOL = np.array(["#1B9E77", "#E7298A", "#a6761d", "#66A61E", "#E6AB02",
                 "#7570B3"])  # metallicities

# Random seed for starlight grid file
RANDSEED = 999
# Reddening law to use in starlight
REDLAW = "CCM"

# Speed of light in km/s
ckms = c.to("km/s").value

# Make plots prettier
rc('font', **{'family': 'serif', 'serif': ['Times New Roman'], 'size': 14})
rc('text', usetex=True)
rc('image', cmap='viridis')


class IFUCube(object):
    """
    Parent to be called by initialisation of child.

    Parameters
    ----------
    cube_hdu : :class:`astropy.io.fits.HDUList`
        A length-3 list of the [primary header, data cube, stddev cube] HDUs
        of the cube to be analysed.
    base_name : str
        The base name to be used for all output.
    RV : float, optional
        The RV value to use for the CCM extinction law.
    sl_dir : None or str, optional
        The directory containing starlight files and bases. The default
        ``None`` will use the `starlight/` subdirectory.
    el_json : None dict or str, optional
        Either a dictonary of emission lines to fit or the filepath to an
        emission lines json file. If ``None`` will use the default located in
        `data/emission_lines.json`. Follow the format of this default to create
        custom line lists.
    """
    def __init__(self, cube_hdu, base_name, RV=3.1, sl_dir=None, el_json=None):

        self.prim_cube = cube_hdu[0]
        self.data_cube = cube_hdu[1]
        self.stddev_cube = cube_hdu[2]

        self.base_name = os.path.abspath(base_name)

        self.RV = RV
        self.data_shape = self.data_cube.data.shape
        self.stddev_shape = self.stddev_cube.shape
        if self.data_shape != self.stddev_shape:
            raise ValueError("data_cube and stddev_cube must have same shape!")

        # Make a wavelength array and store its units
        sl = self.data_cube.header["CRVAL3"]  # start lambda
        rl = self.data_cube.header["CRPIX3"]  # ref pix lambda
        dl = self.data_cube.header["CD3_3"]  # delta lambda
        self.lamb = sl + (np.arange(self.data_shape[0]) - (rl - 1)) * dl
        self.delta_lamb = dl
        self.lamb_unit = u.Unit(self.data_cube.header["CUNIT3"])
        bunit = self.data_cube.header["BUNIT"]
        try:
            self.flux_unit = u.Unit(bunit).to_string("latex")
        except ValueError:
            self.flux_unit = bunit

        if sl_dir is None:
            self.sl_dir = os.path.join(FILEDIR, "starlight")
        elif not os.path.isdir(sl_dir):
            raise AttributeError("{} is not an accesible "
                                 "directory".format(sl_dir))
        else:
            self.sl_dir = sl_dir
        if el_json is None:
            el_filepath = os.path.join(FILEDIR, "data", "emission_lines.json")
            with open(el_filepath) as f:
                self._emission_lines = json.load(f)
        elif isinstance(el_json, dict):
            self._emission_lines = el_json
        else:
            self._emission_lines = json.load(el_json)

        self._z = self.prim_cube.header["IFU_Z"]
        self.nucleus = None
        self.n_cpu = int(min(mp.cpu_count()-1, mp.cpu_count()*0.9))

        # dictionary of {[bin number]: ["starlightoutfile", "spec infile",}
        self.sl_output = {}
        self.results = {}  # Will store (almost) all you could possibly want

    def deredden(self):
        """
        Correct data and stddev cubes for Galactic dust.

        Based on CCM law and the value of the header card `IFU_EBV`.
        After correction, the header card `IFU_EBV` is set to zero.

        """
        ebv = self.prim_cube.header["IFU_EBV"]
        if ebv == 0:
            print("ebv = 0, skipping deredden()")
            return
        print("dereddening with E(B-V) = {:.3f} mag and RV = {}"
              .format(ebv, self.RV))
        Alamb = get_Alamb(self.lamb, ebv, self.RV, self.lamb_unit)[0]
        corr = 10**(0.4 * Alamb)

        # Multiply our flux and stddevs by the correction
        self.data_cube.data *= corr[:, None, None]
        self.stddev_cube.data *= corr[:, None, None]
        # The cube is now dereddened
        self.prim_cube.header["IFU_EBV"] = 0.

    def deredshift(self):
        """
        Correct wavelengths of the data and stddev cubes to restframe.

        After correction, the header card `IFU_Z` is set to zero.
        """
        z = self.prim_cube.header["IFU_Z"]
        if z == 0:
            print("redshift = 0, skipping deredshift()")
            return
        print("deredshifting from z = {}".format(z))
        # Update the lambda array
        self.lamb /= 1. + z
        self.delta_lamb /= 1 + z
        # Update data header
        self.data_cube.header["CRVAL3"] /= 1. + z  # start lambda
        self.data_cube.header["CRPIX3"] /= 1. + z  # ref pix lambda
        self.data_cube.header["CD3_3"] /= 1. + z  # delta lambda
        # Update stddev header
        try:
            self.stddev_cube.header["CRVAL3"] /= 1. + z
            self.stddev_cube.header["CRPIX3"] /= 1. + z
            self.stddev_cube.header["CD3_3"] /= 1. + z
        except KeyError:
            print("warning: couldn't update stddev header wavelength scale")
        # The cube is now restframe
        self.prim_cube.header["IFU_Z"] = 0.

    def mask_regions(self, centres, r):
        """
        Set spaxels in circular regions defined by ``centres`` with radius
        ``r`` to NaNs.

        Parameters
        ----------
        centres : list of lists
            Coordinates of mask centres of the form ((x1,y1),...(xn,yn)).
            Should be given as zero-indexed integer values.
        r : float or array_like
            The radius (radii) of the masks. If given as ``int``, same radius
            used for all masks given in ``centres``. If ``array_like``, this
            will be cycled over for all ``centres``.
        """
        print("masking regions")
        if isinstance(r, (int, float)):
            r = [r]
        elif len(r) != len(centres) and len(r) != 1:
            warnings.warn("Number of radii ({}} not same as centres ({}). "
                          "Radii will be cycled over.", RuntimeWarning)
        for centre, _r in zip(centres, cycle(r)):
            x_cen, y_cen = centre
            y, x = np.ogrid[-y_cen:self.data_shape[1]-y_cen,
                            -x_cen:self.data_shape[2]-x_cen]
            mask = x*x + y*y <= _r*_r
            self.data_cube.data[:, mask] = np.nan

    def set_nucleus(self, xc, yc, usewcs=False, xs=5., ys=5., box_size=10,
                    lamb_low=5500., lamb_upp=5700., plot=True):
        """
        Find the location of the nucleus based on the peak of a gaussian fit.

        Cuts a subarray from the data cube centred at ``xc``, ``yc`` and sums
        the flux for each spaxel in the range ``lamb_low`` to
        ``lamb_upp``. This 2D summation is fitted with a gaussian to determine
        the centre in the image plane. A location can be forced (useful if it
        is outside the FOV, for example) by setting ``box_size`` to 0.  The
        centres can be passed as RA and DEC with ``usewcs``. In this case they
        should be floats or string representations that are parsable by astropy
        `SkyCoord <http://docs.astropy.org/en/stable/coordinates>`_.

        Parameters
        ----------
        xc, yc : float or string, optional
            Best guess for centre of guassian in pixels or RA and DEC
            (see ``usewcs``).
        usewcs : bool, optional
            With ``usewcs``, ``xc`` and ``yc`` will be treated as RA and DEC in
            units of hourangle and degrees, respectively.
        xs, ys : float, optional
            Best guess for sigma of gaussian in pixels.
        box_size : int, optional
            The half-length of the box, centred on ``xc``, ``yc``, to fit with
            a 2D gaussian in pixels. If 0 then set nucleus directly as ``xc``,
            ``yc``.
        lamb_low, lamb_upp : float, optional
            The wavelength limits over which to sum the cube.
        plot : bool, optional
            Make a plot of the data, model and residual).
        """

        if usewcs:
            w = wcs.WCS(self.data_cube.header)
            s = SkyCoord(xc, yc, unit=(u.hourangle, u.deg))
            xc, yc = s.to_pixel(w)
            if np.any(np.isnan([xc, yc])):
                raise ValueError("Couldn't find pixel location for {}"
                                 .format(s.to_string("hmsdms")))
        if isinstance(xc, str) or isinstance(yc, str):
            print("coordinates given as string: if giving RA/DEC, set "
                  "usewcs=True")
            return

        in_x = box_size <= xc <= self.data_shape[2] - box_size
        in_y = box_size <= yc <= self.data_shape[1] - box_size
        if box_size == 0:
            self.nucleus = np.array((round(xc, 3), round(yc, 3)))
            print("set nucleus as {}".format(self.nucleus))
            self.results["nucleus"] = self.nucleus
            return
        elif not in_x or not in_y:
            raise ValueError("box must be fully within the image, use "
                             "box_size=0 to force a location outside "
                             "the FOV.")
            return

        # Find the nearest sampled wavelengths to our limits
        idx_low = np.abs(self.lamb - lamb_low).argmin()
        idx_upp = np.abs(self.lamb - lamb_upp).argmin() + 1

        # Cutout the subarray then sum and normalise it
        cutout = self.data_cube.data[:, yc-box_size:yc+box_size+1,
                                     xc-box_size:xc+box_size+1]
        z = np.sum(cutout[idx_low:idx_upp], axis=0)
        z /= np.max(z)
        y, x = np.mgrid[:z.shape[0], :z.shape[1]]

        # Initial gaussian guess
        g_init = models.Gaussian2D(amplitude=1., x_mean=z.shape[0]/2.,
                                   y_mean=z.shape[1]/2., x_stddev=xs,
                                   y_stddev=ys)
        # Fit the data based on this initial guess
        fit_g = fitting.LevMarLSQFitter()
        g = fit_g(g_init, x, y, z)

        self.nucleus = np.array((round(xc-box_size+g.x_mean, 3),
                                round(yc-box_size+g.y_mean, 3)))
        print("found nucleus as {}".format(self.nucleus))

        self.results["nucleus"] = self.nucleus

        # Plot the data with the best-fit model
        if plot:
            plt.close("all")
            plt.figure(figsize=(9.5, 3))
            for i, vt in enumerate([(z, "datacube"),
                                    (g(x, y), "model"),
                                    (z-g(x, y), "residual")], 1):
                vals, title = vt
                plt.subplot(1, 3, i)
                plt.imshow(vals, origin='lower', interpolation='nearest',
                           vmin=np.min(z), vmax=1.)
                plt.autoscale(False)
                plt.plot(g.x_mean, g.y_mean, "kx", markersize=10)
                plt.gca().get_xaxis().set_visible(False)
                plt.gca().get_yaxis().set_visible(False)
                plt.title(title)
            plt.gcf().text(0.5, 0.05, "nucleus = "+str(self.nucleus),
                           ha="center", va="center")
            plt.savefig(self.base_name+"_nucleus.pdf", bbox_inches="tight")

    def voronoi_bin(self, target_sn, lamb_low=5590., lamb_upp=5680.,
                    cont_lamb_low=None, cont_lamb_upp=None, clobber=False,
                    min_sn=3, weighted=False):
        """
        Apply the voronoi binning algorithm to the data cube.

        The target signal-to-noise (S/N) is given by ``target_sn``. S/N
        calculation is performed in the wavelength window given by ``lamb_low``
        and ``lamb_upp``. Only spaxels that have S/N >= ``min_sn`` will be
        considered for binning. This will write results in a file with the
        suffix `_voronoibins.txt`.

        Emission line signal to noise can be calculated by specifying
        ``cont_lamb_low`` and ``cont_lamb_upp``. In this case the S/N in this
        window is removed from the ``lamb_low`` and ``lamb_upp`` window to give
        the S/N in the emission lines.

        Parameters
        ----------
        target_sn : float, optional
            The target signal-to-noise of the voronoi binning algorithm.
        lamb_low, lamb_upp : float, optional
            The wavelength limits over which to calculate the S/N.
        cont_lamb_low, cont_lamb_upp : float or None optional
            The wavelength limits of the continuum to subtract in calculating
            S/N (defaults to ``None``, ``None`` - i.e. no subtraction)
        clobber: bool, optional
            Whether to overwrite any existing binning results.
        min_sn : float, optional
            The minimum S/N of a spaxel to be considered for the algorithm.
        weighted: bool, optional
            Store the weighted-mean sum of the individual spaxels in each
            bin. If ``false`` just return the arithmetic sum.

        Example
        -------
        To produce bins determined by their signal to noise in H\
        :math:`\\alpha` + [NII]:

        ``>>> IFUCube.voronoi_bin(lamb_low=6540., lamb_upp=6580.,``\
        ``cont_lamb_low=6600, cont_lamb_upp=6640)``

        References
        ----------
        Uses the algorithm of [CC03]_.

        .. [CC03] Cappellari, M. & Copin, Y., "Adaptive spatial binning of
           integral-field spectroscopic data using Voronoi tessellations",
           MNRAS, 2003.

        """
        try:
            self.results["bin"]
        except KeyError:
            pass
        else:
            if not clobber:
                print("bins already exist, use clobber=True to overwrite")
                return

        if self.nucleus is None:
            print("run `set_nucleus` first")
            return

        print("binning spaxels with Voronoi algorithm with "
              "S/N target of {}".format(target_sn))

        if not all((self.lamb[0] < lamb_low < self.lamb[-1],
                    self.lamb[0] < lamb_upp < self.lamb[-1])):
            raise ValueError("S/N window not within wavelength range of cube")
        # Find the nearest sampled wavelengths to our limits
        idx_low = np.abs(self.lamb - lamb_low).argmin()
        idx_upp = np.abs(self.lamb - lamb_upp).argmin() + 1

        # Sum the data and stddev cubes to get signal and noise while
        # catching the warning when spaxels only have nans
        signal = np.nanmean(self.data_cube.data[idx_low:idx_upp, :, :],
                            axis=0)
        noise = np.nanmean(self.stddev_cube.data[idx_low:idx_upp, :, :],
                           axis=0)
        if None not in (cont_lamb_low, cont_lamb_upp):
            cont_idx_low = np.abs(self.lamb - cont_lamb_low).argmin()
            cont_idx_upp = np.abs(self.lamb - cont_lamb_upp).argmin() + 1
            cont = np.nanmean(self.data_cube.data[cont_idx_low:cont_idx_upp,
                                                  :, :], axis=0)
            signal = signal - cont

        # Get the number of nans per spaxel, if its more than 20% of the
        # wavelength window we take this as a crap spectrum and throw it out
        num_nans = np.sum(np.isnan(self.data_cube.data[idx_low:idx_upp, :, :]),
                          axis=0)
        idx_nans = np.where(num_nans > 0.2 * (idx_upp - idx_low))
        signal[idx_nans] = np.nan
        noise[idx_nans] = np.nan

        # Make array to hold x,y coordinates of all spaxels
        xx, yy = np.meshgrid(np.arange(signal.shape[1]),
                             np.arange(signal.shape[0]))
        # Make an array to hold the voronoi input
        vor_input = np.column_stack((xx.ravel(), yy.ravel(),
                                     signal.ravel(), noise.ravel()))
        # Need to clean the input of the bad spectra:
        #  remove those with signal == nan
        vor_input = vor_input[~np.isnan(vor_input[:, 3])]
        #  remove any with negative or zero noise
        vor_input = vor_input[vor_input[:, 3] > 0]
        #  also reduce to only spaxels with S/N >= min_sn
        vor_input = vor_input[(vor_input[:, 2]/vor_input[:, 3]) >= min_sn]

        # Split columns as voronoi wants separate inputs
        x, y, sig, noi = vor_input[:, [0, 1, 2, 3]].T

        # Call the voronoi binning script
        vor_plot = self.base_name + "_bins_voronoi.pdf"
        res = voronoi.voronoi_2d_binning(x, y, sig, noi, targetSN=target_sn,
                                         cvt=True, pixelsize=1, plot=vor_plot,
                                         quiet=False, n_cpu=self.n_cpu)
        bin_num, x_node, y_node, x_bar, y_bar, bin_sn, n_pix, scale = res

        vor_output = np.column_stack([x, y, x_bar[bin_num], y_bar[bin_num],
                                      bin_num])
        # Populate the bin_nums dict
        bin_nums = {}
        for i, bn in enumerate(np.sort(np.unique(bin_num)).astype("int"), 1):
            print("processing bin {}/{}".format(i, len(np.unique(bin_num))),
                  end="\r")
            idx = vor_output[:, 4] == bn
            x, y, x_mean, y_mean, _bn = vor_output[idx].T
            x, y = x.astype("int"), y.astype("int")
            x_mean, y_mean = x_mean[0], y_mean[0]
            nx, ny = self.nucleus - 0.5
            distances = ((x - nx)**2 + (y - ny)**2)**0.5
            spec = self._get_bin_spectrum((x, y), weighted=weighted)
            bin_nums[bn] = {"spax": (x, y),
                            "mean": (x_mean, y_mean),
                            "spec": spec,
                            "dist_min": np.min(distances),
                            "dist_max": np.max(distances),
                            "dist_mean": ((x_mean - nx)**2
                                          + (y_mean - ny)**2)**0.5,
                            "continuum": {"bad": 0},
                            "emission": {"bad": 0}}

        self.results["bin"] = bin_nums
        print()
        print("found {} bins".format(len(bin_nums)))

    def emission_line_bin(self, min_peak_flux, min_frac_flux, max_radius,
                          min_flux, min_npix=8, line_lamb=6562.8, border=3,
                          smooth=0.5, plot=True, clobber=False, weighted=False,
                          **kwargs):
        """
        Apply the HII explorer [SFS]_ binning algorithm to the datacube.

        This method will bin spaxels by attempting to determine distinct
        emission line regions. An emission line map (usually
        H\ :math:`\\alpha`) is created by :func:`~ifuanal.get_line_map` through
        subtraction of a continuum from a narrow band filter centred on the
        emission line. Peaks above ``min_peak_flux`` in the emission line map
        are seeds for bins which are expanded so long as neighbouring spaxels
        are above ``min_frac_flux`` of the peak and within ``max_radius``
        distance. Pixels below ``min_flux`` are excluded from binning. Values
        should be passed in datacube counts. Islands of pixels separated from
        their peak and bins with less than ``min_npix`` pixels are removed

        See :func:`~ifuanal.get_line_map` for more information on the kwargs
        used to define the wavelength windows of the line and continuum and
        their defaults.

        Parameters
        ----------
        min_peak_flux : float
            The minimum flux for a spaxel to be considered as a new bin seed.
        min_frac_flux : float
            The minimum flux of a spaxel, as a fraction of the bin's peak flux,
            to be considered a member of the bin.
        max_radius : float
            The maximum radius allowed for a bin. Any peaks found within
            ``max_radius`` of the image edge will not be used.
        min_flux : float
            The minimum flux of a spaxel for consideration in the binning
            routine.
        min_npix : float, optional
            The minimum number of pixels required for a bin.
        line_lamb : float, optional
            The wavelength of emission line to use (defaults to
            H\ :math:`\\alpha`).
        border : int, optional
            The size of a border in pixels around all nans and the cube edges
            to exclude peaks from. i.e. and peak found within ``border`` of the
            edge of field of view or a masked region will not be used as a bin
            seed.
        smooth : float, optional
            The width of the gaussian filter to use when smoothing the emission
            line map, prior to peak detection, set to zero to skip smoothing.
        plot : bool, optional
            Whether to make a plot of the continuum, filter, line and bin maps.
        clobber : bool, optional
            Whether to overwrite any existing binning results.
        weighted: bool, optional
            Store the weighted-mean sum of the individual spaxels in each
            bin. If ``false`` just return the arithmetic sum.
        filter_width : float, optional
            The filter width to extract around the emission line.
        cont_width : float, optional
            The width of the wavelength window to extract either side of the
            filter to define the continuum.
        cont_pad : float, optional
            The padding to apply between the edge of the filter and the start
            of continuum window.

        References
        ----------
        .. [SFS] S.F. Sanchez, HII_explorer,
           http://www.caha.es/sanchez/HII_explorer/
        """
        try:
            self.results["bin"]
        except KeyError:
            pass
        else:
            if not clobber:
                print("bins already exist, use clobber=True to overwrite")
                return

        if min_flux >= min_peak_flux:
            print("min_peak_flux > min_flux is required")
            return

        if self.nucleus is None:
            print("run `set_nucleus` first")
            return

        print("binning spaxels using HII explorer algorithm around emission "
              "line {}".format(line_lamb))

        r = int(math.ceil(max_radius))

        maps = get_line_map(self.data_cube, self.lamb, line_lamb, **kwargs)
        filt_map, cont_map, line_map = maps

        # Find peaks in the line_map, these can be quite close together
        # as we will grow/merge them later with the HII explorer algorithm
        if smooth > 0:
            gauss = Gaussian2DKernel(stddev=smooth)
            smooth_line_map = convolve(line_map, gauss,
                                       boundary="extend")
        else:
            smooth_line_map = line_map
        line_map_max = ndimage.maximum_filter(smooth_line_map, size=3,
                                              mode="constant")
        # Get the location of the peaks and apply the minimum peak threshhold
        peaks = ((line_map_max == smooth_line_map)
                 & (line_map >= min_peak_flux))
        # Remove all peaks within ``border`` of the cube edge
        peaks[:border, :] = 0
        peaks[-border:, :] = 0
        peaks[:, :border] = 0
        peaks[:, -border:] = 0
        # Get locations in the map of the peaks and sort decending in flux
        sort_idx = line_map[peaks].argsort()[::-1]
        peak_xy = np.argwhere(peaks)[sort_idx]
        # Number of peaks
        n_peaks = len(peak_xy)
        # Keep track of which pixels are allocated to which bin
        bin_map = np.full(peaks.shape, np.nan)
        bin_nums = {}
        bn = 0  # ensure we have contiguous sequence of bin_numbers
        # Starting with the brightest, run the HII explorer algorithm to the
        # map. If a peak is merged with a brighter nearby peak, it is skipped.
        for i, xy in enumerate(peak_xy, 1):
            x, y = map(int, xy)
            # Check if we've already covered this peak
            if ~np.isnan(bin_map[x, y]):
                # bin_peak_xy[i-1] = np.nan
                continue
            print("processing bin seed {}/{}".format(i, n_peaks), end="\r")
            # Check if we're close to nans (i.e. image border or masked
            # region)
            line_cutout = line_map[max(x-border, 0):x+border+1,
                                   max(y-border, 0):y+border+1]
            if np.any(np.isnan(line_cutout)):
                continue
            peak_flux = line_map[x, y]
            thresh_flux = max(peak_flux * min_frac_flux, min_flux)

            # Make cutouts around our peak to check:
            #  the flux of the spaxels
            fluxes = line_map[max(x-r, 0):x+r+1, max(y-r, 0):y+r+1]
            #  if the spaxels have been previously assigned a bin
            allocs = bin_map[max(x-r, 0):x+r+1, max(y-r, 0):y+r+1]
            #  the distance of the spaxels (accounting for if we are close
            #  to the edge of the array)
            x_cen = np.arange(max(-r, -x), min(r+1, line_map.shape[0]-x), 1)
            y_cen = np.arange(max(-r, -y), min(r+1, line_map.shape[1]-y), 1)
            xx, yy = np.meshgrid(y_cen, x_cen)
            dists = np.sqrt(xx * xx + yy * yy)

            # If the spaxel passes these tests it is added to the bin
            bin_cutout = ((dists <= max_radius)
                          & (fluxes >= thresh_flux)
                          & (np.isnan(allocs)))
            # We only want adjoining pixels to be considered so remove
            # separate (not-touching) objects
            struct = ndimage.generate_binary_structure(2, 1)
            objs, nobj = ndimage.label(bin_cutout, structure=struct)
            bin_cutout[np.where(objs != objs[r, r])] = 0
            if np.sum(bin_cutout) >= min_npix:
                # Update the bin_map with this bin number
                bin_map[max(x-r, 0):x+r+1, max(y-r, 0):y+r+1][bin_cutout] = bn
                # Add a bin entry to our dict
                # We swap the x and y to FITS standard in the dict
                xy_spax = np.where(bin_map == bn)[::-1]
                xy_mean = (x, y)[::-1]
                nx, ny = self.nucleus - 0.5
                distances = ((xy_spax[0] - nx)**2
                             + (xy_spax[1] - ny)**2)**0.5
                spec = self._get_bin_spectrum(xy_spax, weighted=weighted)
                bin_nums[bn] = {"spax": xy_spax,
                                "mean": xy_mean,
                                "spec": spec,
                                "dist_min": np.min(distances),
                                "dist_max": np.max(distances),
                                "dist_mean": ((xy_mean[0] - nx)**2
                                              + (xy_mean[1] - ny)**2)**0.5,
                                "continuum": {"bad": 0},
                                "emission": {"bad": 0}}
                bn += 1  # update bin number

        if plot:
            plt.close()
            binfig, ax = plt.subplots(1, 4, sharex=True, sharey=True,
                                      figsize=(16, 4),
                                      subplot_kw={"adjustable": "box-forced"})
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                ax[0].imshow(filt_map, origin="lower",
                             interpolation="none", norm=colors.PowerNorm(0.25))
                ax[0].set_title("Filter")
                ax[0].autoscale(False)
                ax[1].imshow(cont_map, origin="lower",
                             interpolation="none", norm=colors.PowerNorm(0.25))
                ax[1].set_title("Continuum")
                ax[1].autoscale(False)
                ax[2].imshow(line_map, origin="lower",
                             interpolation="none", norm=colors.PowerNorm(0.25))
                ax[2].set_title("Line")
                ax[2].autoscale(False)
                ax[3].imshow(bin_map, origin="lower", interpolation="none",
                             cmap="viridis_r", norm=colors.PowerNorm(0.25))
                ax[3].autoscale(False)
            ax[2].scatter(peak_xy[:, 1], peak_xy[:, 0], s=9, marker="x",
                          c="w", alpha=0.5, lw=0.3)
            ax[3].set_title("Bins")
            binfig.suptitle("Filter centre at {}\\AA".format(line_lamb))
            binfig.savefig(self.base_name+"_bins_el.pdf", bbox_inches="tight")

        self.results["bin"] = bin_nums
        print()
        print("found {} bins".format(len(bin_nums)))

    def nearest_bin(self, min_peak_flux, max_radius, min_flux, min_frac_flux=0,
                    min_npix=8, line_lamb=6562.8, dist_weighted=True,
                    weight_pow=1/3., niter=1, max_filter_size=5, border=3,
                    smooth=0.5, plot=True, clobber=False, weighted=False,
                    **kwargs):
        """
        Create bins based on pixels proximity to peaks, largely derived from
        [LRN]_.

        An emission line map (usually H\ :math:`\\alpha`) is created by
        :func:`~ifuanal.get_line_map` through subtraction of a continuum from a
        narrow band filter centred on the emission line. Peaks above
        ``min_peak_flux`` in the emission line map are seeds for bins. Any peak
        within ``niter`` pixels of a pixel with value less than ``min_flux``
        or a masked (nan) pixel are rejected.  The closest peak to each pixel
        is determined forming the basis of the bins. Those pixels further than
        ``max_radius`` from a peak or below ``min_flux`` are then removed. If
        ``min_frac_flux`` is > 0 then pixels less than the ``min_frac_flux``
        times the peak flux are also removed. Islands of pixels separated from
        the bin underlying their nearest peak, and bins with less than
        ``min_npix`` pixels are removed.

        The optional argument ``weighted`` will apply a weighting to the
        distances of each peak based on the sum of the flux in each unweighted
        bin to the power ``weight_pow``. This means that brighter regions have
        more influence over local pixels than nearby faint peaks. The
        ``max_radius`` argument will still be adhered to.

        See :func:`~ifuanal.get_line_map` for more information on the kwargs
        used to define the wavelength windows of the line and continuum and
        their defaults.

        Parameters
        ----------
        min_peak_flux : float
            The minimum flux for a spaxel to be considered as a new bin seed.
        max_radius : float
            The maximum radius allowed for a pixel to be considered part of a
            bin.
        min_flux : float
            The minimum flux of a spaxel for consideration in the binning
            routine.
        min_frac_flux : float, optional
            The minimum flux of a spaxel, as a fraction of the bin's peak flux,
            to be considered a member of the bin.
        min_npix : float, optional
            The minimum number of pixels required for a bin.
        line_lamb : float, optional
            The wavelength of emission line to use (defaults to
            H\ :math:`\\alpha`).
        dist_weighted : bool, optional
            Whether to weight the distance determination based on the fluxes
            around each peak.
        weight_pow : float, optional
            The power index used in the weighting of distances, requires
            ``dist_weighted = True``
        niter : int, optional
            Peaks within this number of pixels of a masked (nan) or <
            ``min_flux`` pixel will be removed. This can help eliminate noise
            spikes as bin seeds.
        max_filter_size : int, optional
            The square footprint in pixels from which to determine peaks.
            See ``scipy.ndimage.filters.maximum_filter`` documentation.
        border : int, optional
            The size of a border in pixels around all nans and the cube edges
            to exclude peaks from. i.e. and peak found within ``border`` of the
            edge of field of view or a masked region will not be used as a bin
            seed.
        smooth : float, optional
            The width of the gaussian filter to use when smoothing the emission
            line map, prior to peak detection, set to zero to skip smoothing.
        plot : bool, optional
            Whether to make a plot of the continuum, filter, line and bin maps.
        clobber : bool, optional
            Whether to overwrite any existing binning results.
        weighted: bool, optional
            Store the weighted-mean sum of the individual spaxels in each
            bin. If ``false`` just return the arithmetic sum.
        filter_width : float, optional
            The filter width to extract around the emission line.
        cont_width : float, optional
            The width of the wavelength window to extract either side of the
            filter to define the continuum.
        cont_pad : float, optional
            The padding to apply between the edge of the filter and the start
            of continuum window.

        References
        ----------
        .. [LRN] Rousseau-Nepton et al., `NGC628 with SITELLE : I. Imaging
           Spectroscopy of 4285 HII region candidates`, arXiv:1704.05121
        """
        try:
            self.results["bin"]
        except KeyError:
            pass
        else:
            if not clobber:
                print("bins already exist, use clobber=True to overwrite")
                return

        if min_flux >= min_peak_flux:
            print("min_peak_flux > min_flux is required")
            return

        if self.nucleus is None:
            print("run `set_nucleus` first")
            return

        print("binning spaxels using Nearest pixel algorithm around emission "
              "line {}".format(line_lamb))

        print("finding peaks")
        # Get emission line map
        maps = get_line_map(self.data_cube, self.lamb, line_lamb, **kwargs)
        filt_map, cont_map, line_map = maps
        # Gaussian smooth if needed
        if smooth > 0:
            gauss = Gaussian2DKernel(stddev=smooth)
            smooth_line_map = convolve(line_map, gauss,
                                       boundary="extend")
        else:
            smooth_line_map = line_map
        # Find peaks in the line_map
        line_map_max = ndimage.maximum_filter(smooth_line_map,
                                              size=max_filter_size,
                                              mode="constant")
        # Get the location of the peaks and apply the minimum peak threshhold
        peaks = ((line_map_max == smooth_line_map)
                 & (line_map >= min_peak_flux))
        # Remove all peaks within ``border`` of the cube edge
        peaks[:border, :] = 0
        peaks[-border:, :] = 0
        peaks[:, :border] = 0
        peaks[:, -border:] = 0
        # Remove peaks within niter of a pixel < min_flux
        struct = ndimage.generate_binary_structure(2, 1)
        toofaint = ndimage.binary_dilation(line_map < min_flux,
                                           structure=struct,
                                           iterations=niter)
        peaks[toofaint] = False
        sort_idx = line_map[peaks].argsort()[::-1]  # sort descending in flux
        peak_xy = np.argwhere(peaks)[sort_idx]
        n_peaks = len(peak_xy)
        if n_peaks == 0:
            print("no peaks found!")
            return

        print("calculating pixel distances")
        # Distances between pixels and all peaks
        yl, xl = line_map.shape
        x = np.arange(xl)
        y = np.arange(yl)
        yy, xx = np.meshgrid(x, y)
        px_cen = np.column_stack((xx.ravel(), yy.ravel()))
        cdres = cdist(px_cen, peak_xy)
        # Get nearest peaks for each pixel
        nearest = np.argmin(cdres, axis=1)
        # We use zero to mark rejected pixels so the first bin needs to be 1
        near_map = nearest.reshape(line_map.shape).astype(float) + 1
        # Remove nans and <min_flux in original image
        near_map[line_map < min_flux] = 0
        near_map[np.isnan(line_map)] = 0
        # Remove pixels more than max_radius away from any peak
        nearestdist = np.min(cdres, axis=1)
        toofar = (nearestdist > max_radius).reshape(line_map.shape)
        near_map[toofar] = 0
        # Remove pixels not connected to any peak
        labeled_nearmap, nlabels = ndimage.label(near_map)
        good_labels = np.unique(labeled_nearmap[peaks])
        goodlabeled_nearmap = (np.in1d(labeled_nearmap, good_labels)
                               .reshape(labeled_nearmap.shape))
        near_map[~goodlabeled_nearmap] = 0
        # Perform grey closing to remove pixel 'holes' in the bins
        near_map = ndimage.morphology.grey_closing(near_map, size=(3, 3))
        # Mask nans and still abide by distance
        near_map[np.isnan(line_map) | toofar] = 0

        # Weight distances to reform bins if needed
        if dist_weighted:
            print("weighting distances")
            # Calculate fluxes from un-weighted bins
            fluxes = ndimage.sum(line_map, near_map, range(1, n_peaks+1))
            # Recalculate distances, weighting by flux**weight_pow
            cdresw = cdres / fluxes[np.newaxis, :]**weight_pow
            # Add upper limit as max_radius to stop regions 'bleeding'
            dist_mask = cdres > max_radius
            cdresw[dist_mask] = np.inf
            nearestw = np.argmin(cdresw, axis=1)
            near_mapw = nearestw.reshape(line_map.shape).astype(float) + 1
            near_mapw[near_map == 0] = 0  # apply all bin rejection as above
            near_map = near_mapw
        # Remove pixels not connected to *their* peak or below min_frac_flux
        # of the peak
        for i, xy in enumerate(peak_xy, 1):
            _nrmap = np.copy(near_map)
            _nrmap[_nrmap != i] = 0
            _lblmap, _nlbl = ndimage.label(_nrmap)
            _goodlbl = np.unique(_lblmap[peaks])
            _goodlblmap = np.in1d(_lblmap, _goodlbl).reshape(_nrmap.shape)
            near_map[(~_goodlblmap) & (_nrmap > 0)] = 0
            peak_val = line_map[xy[0], xy[1]]
            near_map[(line_map < min_frac_flux * peak_val) & (_nrmap > 0)] = 0
        # Change back to zero-index labels and nans for bad pixels
        near_map[near_map == 0] = np.nan
        near_map -= 1

        # Fill bin_nums dict for storage in results dict
        bin_nums = {}
        for bn, xy in enumerate(peak_xy):
            x, y = map(int, xy)
            print("processing bin {}/{}".format(bn+1, n_peaks), end="\r")
            sys.stdout.flush()
            xy_spax = np.where(near_map == bn)[::-1]
            xy_mean = (x, y)[::-1]
            nx, ny = self.nucleus - 0.5
            distances = ((xy_spax[0] - nx)**2
                         + (xy_spax[1] - ny)**2)**0.5
            spec = self._get_bin_spectrum(xy_spax, weighted=weighted)
            bin_nums[bn] = {"spax": xy_spax,
                            "mean": xy_mean,
                            "spec": spec,
                            "dist_min": np.min(distances),
                            "dist_max": np.max(distances),
                            "dist_mean": ((xy_mean[0] - nx)**2
                                          + (xy_mean[1] - ny)**2)**0.5,
                            "continuum": {"bad": 0},
                            "emission": {"bad": 0}}
        if plot:
            plt.close()
            binfig, ax = plt.subplots(1, 4, sharex=True, sharey=True,
                                      figsize=(16, 4),
                                      subplot_kw={"adjustable": "box-forced"})
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                ax[0].imshow(filt_map, origin="lower",
                             interpolation="none", norm=colors.PowerNorm(0.25))
                ax[0].set_title("Filter")
                ax[0].autoscale(False)
                ax[1].imshow(cont_map, origin="lower",
                             interpolation="none", norm=colors.PowerNorm(0.25))
                ax[1].set_title("Continuum")
                ax[1].autoscale(False)
                ax[2].imshow(line_map, origin="lower",
                             interpolation="none", norm=colors.PowerNorm(0.25))
                ax[2].set_title("Line")
                ax[2].autoscale(False)
                ax[3].imshow(near_map, origin="lower", interpolation="none",
                             cmap="viridis_r", norm=colors.PowerNorm(0.25))
                ax[3].autoscale(False)
            ax[2].scatter(peak_xy[:, 1], peak_xy[:, 0], s=9, marker="x",
                          c="w", alpha=0.5, lw=0.3)
            ax[3].set_title("Bins")
            binfig.suptitle("Filter centre at {}\\AA".format(line_lamb))
            binfig.savefig(self.base_name+"_bins_nr.pdf", bbox_inches="tight")

        self.results["bin"] = bin_nums
        print()
        print("found {} bins".format(len(bin_nums)))

    def add_custom_bin(self, centre, r=None, xy=None, weighted=False):
        """
        Create a custom bin to analyse. Either specify ``centre`` and
        ``r`` to create a circular bin, or give list of x and y coordinates as
        ``xy`` in order to select specific spaxels. For the ``xy`` specified
        bins, a ``centre`` is still required and must be provided.

        Custom bins have negative values (beginning at -1) in all output etc.

        Parameters
        ----------
        centre : array_like
            Length 2 array giving the x,y centre of the bin
        r : int, optional
            The radius of the bin in pixels.
        xy : array_like, optional
            Length 2 array of x and y pixels to use as the bin
        weighted : bool, optional
            Whether to return the weighted-mean sum or arthimetic sum spectrum.

        Example
        -------
        Creating an SDSS-like  3 arcsec fibre centred on the host:

        1. set the nucleus with :meth:`set_nucleus()` (it must be in the FOV)
        2. determine the ``plate_scale`` of the cube (arcsec/pixel)
        3. ``>>> cube.add_custom_bin(cube.nucleus, 3/plate_scale)``
        """
        if self.nucleus is None:
            print("run `set_nucleus` first")
            return

        x_cen, y_cen = centre
        if r is not None:
            y, x = np.ogrid[-y_cen:self.data_shape[1]-y_cen,
                            -x_cen:self.data_shape[2]-x_cen]
            xy_spax = np.where(x*x + y*y <= r*r)[::-1]
        elif xy is not None:
            xy_spax = (np.atleast_1d(xy[0]), np.atleast_1d(xy[1]))
        else:
            print("either specify `r`, or `xy`")

        bn = -1
        try:
            existing_bins = self._get_bin_nums("all")
        except KeyError:
            # we have no bins created yet
            self.results["bin"] = {}
        else:
            while True:
                if bn not in existing_bins:
                    break
                bn -= 1

        nx, ny = self.nucleus - 0.5
        distances = ((xy_spax[0] - nx)**2 + (xy_spax[1] - ny)**2)**0.5
        spec = self._get_bin_spectrum(xy_spax, weighted=weighted)
        self.results["bin"][bn] = {"spax": xy_spax,
                                   "mean": (x_cen, y_cen),
                                   "spec": spec,
                                   "dist_min": np.min(distances),
                                   "dist_max": np.max(distances),
                                   "dist_mean": ((x_cen - nx)**2
                                                 + (y_cen - ny)**2)**0.5,
                                   "continuum": {"bad": 0},
                                   "emission": {"bad": 0}}
        print("added custom bin {} to the list".format(bn))

    def get_loc_bin(self, loc):
        """
        Return the bin number for a given location.

        Parameters
        ----------
        loc : array_like
            Length 2 array giving the (``x``, ``y``) pixel coordinates to get
            the bin number of.

        Returns
        -------
        bin_num : int
            The bin number that contains ``loc``. Returns None if not in any
            bin.
        """
        x, y = map(int, loc)
        for bn, d in self.results["bin"].items():
            bn_x, bn_y = d["spax"]
            if x in bn_x and y in bn_y:
                return bn
        print("didn't find a bin at location {}, {}".format(x, y))
        return None

    def _get_bin_nums(self, bins="nobad", custom=False):
        """
        Get the bin numbers and, optionally, the number of custom bins

        Parameters
        ----------
        bins : str, optional
            Should be 1) "all", 2) "nocontbad", 3) "nobad".  1) returns all
            bins present in results["bin"]. 2) returns only those with good
            continuum fitting. 3) returns only those with good emission and
            continuum fitting.
        custom : bool, optional
            Whether to return the number of custom bins

        Returns
        -------
        bin_nums : list or tuple
            If custom returns tuple of (bin numbers, number custom bins) else
            returns list of bin numbers
        """
        if bins not in ("all", "nocontbad", "nobad"):
            print("bins must be one of ('all', 'nocontbad', 'nobad')")
            return
        bin_nums = np.sort(self.results["bin"].keys())
        if bins in ("nocontbad", "nobad"):
            idx_remove = []
            for i, bn in enumerate(bin_nums):
                contbad = self.results["bin"][bn]["continuum"]["bad"]
                emisbad = self.results["bin"][bn]["emission"]["bad"] * 2
                bad_value = contbad + emisbad
                if bad_value != 0 and bins == "nobad":
                    idx_remove.append(i)
                elif bad_value not in (0, 2) and bins == "nocontbad":
                    idx_remove.append(i)
            bin_nums = np.delete(bin_nums, idx_remove)
        if custom:
            return (bin_nums, np.sum(bin_nums < 0))
        return bin_nums

    def _get_single_spectrum(self, x, y):
        """
        Return the spectrum for spaxel at location ``x``, ``y``.

        Returned array is shape (N,4), where N is number of wavelength
        samples of cube and columns are wavelength, flux, flux_err, flag.
        Flag = 2 for negative or nan fluxes (dummy values for flux and flux_err
        are inserted for these wavelengths).
        """
        if not all((len(x) == 1, len(y) == 1)):
            raise AttributeError("``x`` and ``y`` should be length-1")
        x = int(x)
        y = int(y)
        spec = np.empty((self.lamb.size, 4))
        spec[:, 0] = self.lamb

        spec[:, 1] = self.data_cube.data[:, y, x]  # numpy axes switch
        spec[:, 2] = self.stddev_cube.data[:, y, x]
        # STARLIGHT ignores flags >=2
        bad_idx = np.isnan(spec[:, 1]) | np.isnan(spec[:, 2])
        spec[:, 3] = (bad_idx).astype("int") * 2

        return spec

    def _get_multi_spectrum(self, x, y, weighted=False):
        """
        Return the summed spectrum over multiple spaxels at locations ``x``,
        ``y``.

        Similar to ``_get_single_spectrum`` except ``x`` and ``y`` are arrays.
        The single spectra given by these locations are combined using the
        weighted-mean or arithmetic sum of the fluxes depending on `weighted``.
        Returns array of same form as ``_get_single_spectrum``.

        """
        x = np.ravel(np.asarray(x))
        y = np.ravel(np.asarray(y))
        if x.size != y.size:
            raise AttributeError("``x`` and ``y`` should be same size")
        # Use weighted arthimetic mean and variance of weighted mean
        # arrays to hold all the flux and flux_stddev values in spaxels
        spaxels_flux = self.data_cube.data[:, y, x]  # numpy axes switch
        spaxels_stddev = self.stddev_cube.data[:, y, x]

        # Find any wavelengths where we have >75% spaxels as nans
        # and flag that wavelength as bad for starlight
        bad_idx = np.isnan(spaxels_flux) | np.isnan(spaxels_stddev)
        num_bad = np.sum(bad_idx, axis=1)
        bad_lamb = num_bad > 0.75 * x.size

        # Array to hold final weighted-mean spectrum - same format
        # as _get_single_spectrum()
        spec = np.empty((self.lamb.size, 4))
        spec[:, 0] = self.lamb
        # Use masked arrays to cover the nans while preserving shape
        spaxels_flux_ma = np.ma.masked_array(spaxels_flux, bad_idx)
        spaxels_stddev_ma = np.ma.masked_array(spaxels_stddev, bad_idx)
        if weighted:
            # Calculate the weighted mean and uncertainty
            w = 1/spaxels_stddev_ma**2
            spec[:, 1] = np.ma.average(spaxels_flux_ma, weights=w, axis=1) * x.size
            spec[:, 2] = 1/np.sum(w, axis=1)**0.5 * x.size
        else:
            spec[:, 1] = np.sum(spaxels_flux_ma, axis=1)
            spec[:, 2] = np.sqrt(np.sum(spaxels_stddev_ma**2, axis=1))
        # STARLIGHT ignores flags >=2
        spec[:, 3] = bad_lamb.astype("int") * 2

        return spec

    def _get_bin_spectrum(self, xy, weighted=False):
        """
        Return the spectrum for coordinates xy.

        Calls ``_get_single_spectrum`` or ``_get_multi_spectrum`` as
        appropriate (i.e. if number of spaxels in bin is 1 or >1).

        Parameters
        ----------
        xy : Nx2 array
            The ([x0, x1.. xN], [y0, y1, yN]) pixel coordinates of the spectrum
            to return.
        weighted : bool, optional
            If ``True`` return the weight-mean sum spectrum, otherwise simply the
            arithmetic sum.

        Returns
        -------
        spec : 4xM array
            see :meth:`ifuanal.IFUCube._get_single_spectrum` and
            :meth:`ifuanal.IFUCube._get_multi_spectrum`
        """
        n_spax = len(xy[0])
        if len(xy[1]) != n_spax:
            print("must have same number of x and y coordinates")
            return None
        if n_spax == 0:
            print("no spaxels found!")
            return None
        elif n_spax == 1:
            # yield the single spaxel's spectrum
            spec = self._get_single_spectrum(xy[0], xy[1])
        else:
            spec = self._get_multi_spectrum(xy[0], xy[1], weighted=weighted)
        return spec

    def run_starlight(self, bin_num=None, base_name="bc03", lamb_low=5590.,
                      lamb_upp=5680., tmp_dir=None, clobber=False):
        """
        Run the starlight fitting of spectra in the cube.

        Creates a unique directory with prefix `sl_` in which to store
        starlight output. By default this is in the system's `$TMP`, but can be
        manually set with ``tmp_dir``. starlight doe not like long file paths:
        ensure the directory file path is short (<40 chars or so).

        Parameters
        ----------
        bin_num : (None, int or array_like), optional
            The bin numbers to fit with starlight (defaults to None, this will
            fit all bins).
        base_name : str
            The name of the bases to use for starlight fitting. There must
            exist a valid starlight bases file (see starlight docs) named
            "starlight/``base_name``.base" (defaults to "bc03").
        lamb_low, lamb_upp : float, optional
            The wavelength limits over which starlight will calculate the S/N.
            This isn't actually used since we have a stddev cube but the limits
            must be within the wavelength range of the cube anyway.
        tmp_dir : str, optional
            The directory in which to store the starlight output. A unique
            directory with prefix `sl_` will be created in ``tmp_dir``.
        clobber : bool, optional
            Whether to overwrite pre-existing results.
        """
        print("running starlight fitting")

        if bin_num is None:
            bin_num = self._get_bin_nums("all")
        elif isinstance(bin_num, (int, float)):
            bin_num = [bin_num]
        if not np.all(np.in1d(bin_num, self._get_bin_nums("all"))):
            raise ValueError("one or more bin numbers given do not exist")

        if not clobber:
            for bn in bin_num:
                try:
                    self.results["bin"][bn]["continuum"]["sl_output"]
                except KeyError:
                    pass
                else:
                    print("previous results found. use clobber=True to "
                          "overwrite.")
                    return

        print("fitting {} bins...".format(len(bin_num)))
        # We need to copy everything to a temporary folder since
        # starlight has filepath character length issues, plus
        # then it doesn't clog up cwd.
        # Make a safe temporary directory to hold everything
        sl_tmp_dir = tempfile.mkdtemp(prefix="sl_", dir=tmp_dir)
        sl_tmp_dir = os.path.join(os.path.abspath(sl_tmp_dir), "")
        if len(sl_tmp_dir) > 40:
            warnings.warn("Temporary directory file path is quite long. "
                          "starlight may have problems: {}"
                          .format(sl_tmp_dir), RuntimeWarning)
        print("STARLIGHT tmp directory for this run is {}"
              .format(sl_tmp_dir))
        # Because shutil.copytree requires `dst` doesn't exist:
        shutil.rmtree(sl_tmp_dir)
        shutil.copytree(self.sl_dir, sl_tmp_dir)
        d = os.path.join(sl_tmp_dir, base_name)
        basefiles = [os.path.join(d, f) for f in os.listdir(d)]
        Nbasefiles = len(basefiles)
        for i, basefile in enumerate(basefiles, 1):
            print("resampling base files {:>4}/{:4}".format(i, Nbasefiles),
                  end="\r")
            resample_base(basefile, self.lamb, self.delta_lamb)
        print()
        # Compute the spectra and write the spectrum to a temp file
        spec_files = []
        for i, bn in enumerate(bin_num, 1):
            print("writing spec files {}/{}".format(i, len(bin_num)),
                  end="\r")
            with tempfile.NamedTemporaryFile(prefix="spec_", dir=sl_tmp_dir,
                                             delete=False) as spec_file:
                np.savetxt(spec_file, self.results["bin"][bn]["spec"],
                           fmt="%14.8f")
                spec_files.append(spec_file.name)
        print()
        # multiprocessing params for starlight pool of workers
        n_cpu = min(self.n_cpu, len(bin_num))
        p = mp.Pool(n_cpu)
        bin_out_files = p.map(fit_starlight,
                              zip(bin_num,
                                  spec_files,
                                  repeat(self.lamb),
                                  repeat(self.delta_lamb),
                                  repeat(sl_tmp_dir),
                                  repeat(base_name),
                                  repeat(lamb_low),
                                  repeat(lamb_upp)))
        p.close()
        p.join()
        print("STARLIGHT fits of {} bins complete".format(len(bin_num)))

        for bn, outfile in bin_out_files:
            self.results["bin"][bn]["continuum"]["sl_output"] = outfile

        self._parse_continuum(bin_num)

    def _parse_continuum(self, bin_num=None):
        """
        Parses the output from STARLIGHT and populates the results dict.

        Parameters
        ----------
        bin_num : (None, int or array_like), optional
            The bin numbers to parse STARLIGHT results for (defaults to None,
            this will parse all bins).
        """
        print("parsing results")
        if bin_num is None:
            bin_num = self._get_bin_nums("all")
        elif isinstance(bin_num, (int, float)):
            bin_num = [bin_num]
        if not np.all(np.in1d(bin_num, self._get_bin_nums("all"))):
            print("one or more bin numbers given do not exist")
            return

        n_bins = len(bin_num)
        for i, bn in enumerate(bin_num, 1):
            print("parsing starlight output {:5}/{:5}".format(i, n_bins),
                  end="\r")
            bin_res = self.results["bin"][bn]["continuum"]
            if bin_res["sl_output"] is not None:
                try:
                    open(bin_res["sl_output"])
                except IOError:
                    print("couldn't open {} for bin {}"
                          .format(bin_res["sl_output"], bn))
                    continue
                sl_results = parse_starlight(bin_res["sl_output"])
                if sl_results is None:  # i.e. failed to get Nl_Obs
                    print("failed to parse {} for bin {}"
                          .format(bin_res["sl_output"], bn))
                    bin_res["bad"] = 1
                    continue
                sl_results["ebv_star"] = sl_results["AV_min"]/self.RV
                bin_res["bad"] = 0
                bin_res.update(sl_results)
            else:
                print("no sl_output found for bin {:<20}".format(bn))
                bin_res["bad"] = 1
        print()

    def run_emission_lines(self, bin_num=None, vd_init=[10, 40, 70, 100],
                           v0_init=[-300, -200, -100, 0, 100, 200, 300],
                           amp_init=[0.1, 1, 10], stddev_bounds=[5, 120],
                           offset_bounds=[-500, 500], weights=True,
                           filtwidth=None, clobber=False):
        """
        Fit emission lines in the continuum-subtracted spectra with gaussians.

        A series of gaussians will be fit to the emission lines.  The fitting
        is performed in a brute force manner over a range of initial guesses to
        attempt to find the global minimum. As such, increasing the length of
        the initial guesses for the width, mean and amplitude (``vd_init``,
        ``v0_init`` and ``amp_init`` respectively) will dramatically increase
        computing time.

        Parameters
        ----------
        bin_num : (None, int or array_like), optional
            The bin numbers to fit (defaults to None, this will
            fit all bins where continuum fitting was sucessful)
        vd_init : list, optional
            The intial guess(es) for the velocity dispersion of the lines in
            km/s. If None will use initial guesses of 10, 40, 70 and 100 km/s.
        v0_init : list or "vstar", optional
            The intial guess(es) for the velocity offset of the lines in
            km/s. The special case "vstar" will use the velocty offset
            determined by STARLIGHT and construct 6 guesses over the range
            +/-200 km/s of that value, but limited to within ``offset_bounds``
            of zero offset.
        amp_init : list, optional
            The initial guess(es) for the amplitudes of the lines in units of
            ``fobs_norm`` (see starlight manual).
        stddev_bounds : list, optional
            Bounding limits for the sigma widths of the lines in km/s.
        offset_bounds : list, optional
            Bounding limits for the offset of the mean of the lines in km/s.
        weights : bool, optional
            Whether to include weights (1/stddev) in the fitting and chi2
            determination.
        filtwidth : float, optional
            The size in spectral wavelength units of a median filter kernal to
            use to fit the residuals after continuum subtraction.  Should be
            large enough to not fit to emission lines (>50 or so).  If ``None``
            (default) then this step is skipped. The residual function fit and
            the median filter width as the number of spectral elements (i.e
            filtwidth/delta_lamb) are stored in the bin results as ``resid_fn``
            and ``filtwidth``, respectively.
        clobber : bool, optional
            Whether to overwrite pre-existing results.
        """
        if bin_num is None:
            bin_num = self._get_bin_nums("nocontbad")
            if len(bin_num) == 0:
                print("found no bins with good continuum fitting")
                return
        elif isinstance(bin_num, (int, float)):
            bin_num = [bin_num]
        if not np.all(np.in1d(bin_num, self._get_bin_nums("nocontbad"))):
            print("one or more bins have bad continuum fitting or don't exist")
            return

        if not clobber:
            for bn in bin_num:
                try:
                    self.results["bin"][bn]["emission"]["chi2dof"]
                except KeyError:
                    pass
                else:
                    print("previous results found. use clobber=True to "
                          "overwrite.")
                    return

        print("fitting emission lines to {} bins...".format(len(bin_num)))
        # convert filtwidth from angstroms to N spectral elements
        if filtwidth:
            filtwidth /= self.delta_lamb
        # multiprocessing params for emission line fitting pool of workers
        n_cpu = min(self.n_cpu, len(bin_num))
        p = mp.Pool(n_cpu)
        emline_results = p.map(fit_emission_lines,
                               zip(bin_num,
                                   [self.results["bin"][bn] for bn in bin_num],
                                   repeat(self._emission_lines),
                                   repeat(vd_init),
                                   repeat(v0_init),
                                   repeat(amp_init),
                                   repeat(stddev_bounds),
                                   repeat(offset_bounds),
                                   repeat(weights),
                                   repeat(filtwidth)))
        p.close()
        p.join()
        print()
        print("emission line fitting complete")
        for bn, emline_res in emline_results:
            self.results["bin"][bn]["emission"] = emline_res

        self._parse_emission(bin_num)

    def _parse_emission(self, bin_num=None, cont_order=2, snrlimit=3):
        """
        Calculate useful quantities from the fitted emission line model.

        The values are calculated for each bin in ``bin_num`` and are based
        on the emission line models produced by run_emission_lines(). Fluxes
        are corrected for the balmer decrement estimate of reddening. Where
        this is not possible (i.e. SNR < 3 for Hbeta), the value of ``ebv``
        is set to ``nan`` and no correction applied.

        Parameters
        ----------
        bin_num : (None, int or array_like), optional
            The bin numbers to get parse emission results for (defaults to
            None, this will fit all bins where continuum fitting was sucessful)
        cont_order : int, optional
            The order of the polynomial to fit to the continuum model in order
            to determine the continuum level for equivalent width measurements.
        snrlimit: float, optional
            Metallicity values will only be calculated where all dependant
            lines are detected at a SNR above this value.
        """

        if bin_num is None:
            bin_num = self._get_bin_nums("nocontbad")
        elif isinstance(bin_num, (int, float)):
            bin_num = [bin_num]
        if not np.all(np.in1d(bin_num, self._get_bin_nums("nocontbad"))):
            print("one or more bins have bad continuum fitting or don't exist")

        n_bins = len(bin_num)
        for i, bn in enumerate(bin_num, 1):
            print("parsing emission model {:>5}/{:5}".format(i, n_bins),
                  end="\r")
            bin_res = self.results["bin"][bn]["emission"]
            bin_res_c = self.results["bin"][bn]["continuum"]
            try:
                emlines = bin_res["lines"]
            except KeyError:
                print("skipping bin {}, no emission line results found"
                      .format(bn))
                continue

            # Determine fluxes, EW, offsets, FWHM etc. for each line
            for line, fit in emlines.items():
                amp, mean, stddev = fit["fit_params"]
                amp_sig, mean_sig, stddev_sig = fit["fit_uncerts"]

                # Get the continuum level with a `cont_order` polynomial
                # fitted over a +/- 100AA window about the line and sample
                # at the line's mean and correct for normalisation
                low_idx = max(np.abs(self.lamb - (mean - 100)).argmin(), 0)
                upp_idx = min(np.abs(self.lamb - (mean + 100)).argmin(),
                              len(self.lamb))
                # Remove pixels that were masked/clipped in continuum fitting
                cont_px = bin_res_c["sl_spec"][low_idx:upp_idx, 3] > 0
                if ((np.sum(cont_px[:len(cont_px)/2]) < 10 or
                     np.sum(cont_px[len(cont_px)/2:]) < 10)):
                    # if we don't have enough good continuum fit pixels in the
                    # window then just take everything outside the emission
                    # line mask
                    cont_px = bin_res_c["sl_spec"][low_idx:upp_idx, 3] != 0
                cont_wl = bin_res_c["sl_spec"][low_idx:upp_idx, 0][cont_px]
                cont_fl = bin_res_c["sl_spec"][low_idx:upp_idx, 1][cont_px]
                pcoeff = np.polyfit(cont_wl, cont_fl, cont_order)
                cont_poly = np.poly1d(pcoeff)
                cont = (cont_poly(mean) * bin_res_c["fobs_norm"])
                # Estimate the continuum level uncertainty as the rms of the
                # residuals about this polynomial
                cont_resids = cont_poly(cont_wl) - cont_fl
                cont_uncert = (np.average(cont_resids**2)**0.5
                               * bin_res_c["fobs_norm"])

                # Flux and EW and their uncertainty components
                flux = (2*math.pi)**0.5 * stddev * amp
                ew = flux/cont
                flux_statuncert = ((2*math.pi) * ((stddev_sig**2 * amp**2)
                                            + (stddev**2 * amp_sig**2)))**0.5
                N = 6 * stddev
                flux_photonuncert = (cont_uncert * N**0.5
                                    * (1 + ew/(N * self.delta_lamb))**0.5)
                flux_uncert = (flux_statuncert**2 + flux_photonuncert**2)**0.5
                ew_uncert = ew * ((flux_uncert/flux)**2
                                  + (cont_uncert/cont)**2)**0.5

                emlines[line]["cont"] = [cont, cont_uncert]
                emlines[line]["flux"] = np.array([flux, flux_uncert])
                emlines[line]["ew"] = np.array([ew, ew_uncert])
                emlines[line]["snr"] = flux/flux_uncert

                def to_fwhm(x): return x * 2 * (2*math.log(2))**0.5
                rest = emlines[line]["rest_lambda"]
                emlines[line]["fwhm"] = [to_fwhm(stddev) * ckms / rest,
                                         to_fwhm(stddev_sig) * ckms / rest]
                emlines[line]["offset"] = [(mean - rest) * ckms / rest,
                                           mean_sig * ckms / rest]
                emlines[line]["mean"] = [mean, mean_sig]

            # Calculate E(B-V)_gas based on balmer decrement
            ebv = np.nan
            try:
                Hb_snr = emlines["Hbeta_4861"]["snr"]
                Ha_snr = emlines["Halpha_6563"]["snr"]
            except KeyError:
                pass
            else:
                if (Hb_snr > 3 and Ha_snr > 3):
                    flux_Hb = emlines["Hbeta_4861"]["flux"][0]
                    flux_Ha = emlines["Halpha_6563"]["flux"][0]
                    lamb_Hb_Ha = [emlines["Hbeta_4861"]["rest_lambda"],
                                  emlines["Halpha_6563"]["rest_lambda"]]
                    k_Hb, k_Ha = get_Alamb(lamb_Hb_Ha, 0.0, self.RV)[1]
                    ebv = (2.5/(k_Hb - k_Ha)) * np.log10((flux_Ha/flux_Hb)/2.86)
                    # Correct fluxes for E(B-V)_gas
                    # FIXME uncertainty in E(B-V)_gas not propagated onto new
                    # dereddened fluxes (errors are correlated)
                    for line in emlines:
                        Alamb = get_Alamb(emlines[line]["mean"][0],
                                          ebv, self.RV)[0]
                        corr = 10**(0.4 * Alamb)
                        emlines[line]["flux"] *= corr
                        emlines[line]["ew"] *= corr
            bin_res["ebv_gas"] = ebv

            # Determine metallicities where possible using SNR>3 lines only
            # Make a short of [flux, flux uncert] to reduce clutter and set
            # initially to np.nan to stop KeyError if line is not present
            el = {"Hbeta_4861": (np.nan, np.nan),
                  "[OIII]_5007": (np.nan, np.nan),
                  "Halpha_6563": (np.nan, np.nan),
                  "[NII]_6583": (np.nan, np.nan),
                  "[SII]_6716": (np.nan, np.nan),
                  "[SII]_6731": (np.nan, np.nan),
                  }
            for line, d in emlines.items():
                if d["snr"] > snrlimit:
                    el[line] = d["flux"]
            # N2
            NII, NII_uncert = el["[NII]_6583"]
            Ha, Ha_uncert = el["Halpha_6563"]
            N2f = NII/Ha  # in flux
            N2f_uncert = N2f * ((NII_uncert/NII)**2 + (Ha_uncert/Ha)**2)**0.5
            # O3N2
            OIII, OIII_uncert = el["[OIII]_5007"]
            Hb, Hb_uncert = el["Hbeta_4861"]
            O3f = OIII/Hb
            O3f_uncert = O3f * ((OIII_uncert/OIII)**2 + (Hb_uncert/Hb)**2)**0.5
            O3N2f = O3f/N2f
            O3N2f_uncert = O3N2f * ((O3f_uncert/O3f)**2
                                    + (N2f_uncert/N2f)**2)**0.5
            # S2
            SII = el["[SII]_6716"][0] + el["[SII]_6731"][0]
            SII_uncert = (el["[SII]_6716"][1]**2 + el["[SII]_6731"][1]**2)**0.5
            S2f = NII/SII
            S2f_uncert = S2f * ((NII_uncert/NII)**2 + (SII_uncert/SII)**2)**0.5
            # Log of values
            N2 = np.log10(N2f)
            N2_uncert = 0.434 * (N2f_uncert/N2f)  # uncertainty on log(NII/Ha)
            O3N2 = np.log10(O3f/N2f)
            O3N2_uncert = 0.434 * (O3N2f_uncert/O3N2f)
            S2 = np.log10(S2f)
            S2_uncert = 0.434 * (S2f_uncert/S2f)
            y = S2 + 0.264 * N2
            y_uncert = (S2_uncert**2 + (0.264*N2_uncert)**2)**0.5
            # Metallicity indicators and uncerts
            PP04_N2 = [8.90 + 0.57 * N2, 0.57 * N2_uncert]
            PP04_O3N2 = [8.73 - 0.32 * O3N2, 0.32 * O3N2_uncert]
            M13_N2 = [8.743 + 0.462 * N2, 0.462 * N2_uncert]
            M13_O3N2 = [8.533 - 0.214 * O3N2, 0.214 * O3N2_uncert]
            D16 = [8.77 + y, y_uncert]
            # Assign to dictionary
            bin_res["metallicity"] = {}
            bin_res["metallicity"]["PP04_N2"] = PP04_N2
            bin_res["metallicity"]["PP04_O3N2"] = PP04_O3N2
            bin_res["metallicity"]["M13_N2"] = M13_N2
            bin_res["metallicity"]["M13_O3N2"] = M13_O3N2
            bin_res["metallicity"]["D16"] = D16

    def make_2dfits(self, keys, suffix, idx=0, uncert_idx=None, clobber=False):
        """
        Create a 2D fits image of bin values.

        The choice of value to make a 2D map of is specified by a list of
        ``keys`` in order to navigate each bin's results dictionary. See
        tutorial documentation for layout of ``results`` dict. The ``idx``
        index of the value is used, optionally a second extension can be saved
        which will use the ``uncert_idx`` index of the value (generally the
        uncertainty).

        Parameters
        ----------
        keys : list
            A list of keys in the ``results`` dict to navigate to the desired
            bin value.
        suffix : str
            The suffix to add to the cube's file name when saving the fits
            file. Saved as `[cube name]_2D_`\ ``suffix``\ `.fits`
        idx : int, optional
            If the value is an array, this specifies the index of the array to
            plot.
        uncert_idx : None or int, optional
            If not ``None``, this index of the array will be saved as a second
            extension in the fits file (e.g. to save the uncertainty associated
            with each value).
        clobber : bool, optional
            Whether to overwrite an existing file.

        Example
        -------
        Create a 2D fits image of the H\ :math:`\\alpha` flux in each bin

        ``>>> IFUCube.results_to_2dfits(["emission", "lines", "Halpha_6563",
        "flux"], "Halphaflux")``

        As above but including the uncertainty in a second extension

        ``>>> IFUCube.results_to_2dfits(["emission", "lines", "Halpha_6563",``\
        ``"flux"], "Halphafluxuncert", uncert_idx=1)``

        Create a D16 metallicity map with uncertainty in second extension

        ``>>> IFUCube.results_to_2dfits(["emission", "metallicity", "D16"],
        "D16", uncert_idx=1)``
        """
        outfile = "{}_2D_{}.fits".format(self.base_name, suffix)
        if not clobber and os.path.isfile(outfile):
                print("{} exists and clobber is false".format(outfile))
                return

        twodmap = np.full(self.data_cube.shape[1:], np.nan)
        if uncert_idx is not None:
            twodmap_uncert = np.full(self.data_cube.shape[1:], np.nan)
        if "emission" in keys:
            bin_nums = self._get_bin_nums("nobad")
        elif "continuum" in keys:
            bin_nums = self._get_bin_nums("nocontbad")
        else:
            bin_nums = self._get_bin_nums("all")

        try:
            reduce(operator.getitem, keys, self.results["bin"][bin_nums[0]])
        except KeyError:
            print("Couldn't find values in results dict given keys: ",
                  keys)
            return

        for bn in bin_nums:
            bin_res = self.results["bin"][bn]
            try:
                dict_val = reduce(operator.getitem, keys, bin_res)
            except KeyError:
                print("{} ".format(bn), end="")
            else:
                try:
                    val = dict_val[idx]
                except (IndexError, TypeError):
                    val = dict_val
                twodmap[bin_res["spax"][::-1]] = val
                if uncert_idx is not None:
                    try:
                        val_uncert = dict_val[uncert_idx]
                    except (IndexError, TypeError):
                        print("{}(u) ".format(bn), end="")
                    else:
                        twodmap_uncert[bin_res["spax"][::-1]] = val_uncert
        print()
        hdulist = fits.HDUList()
        wcs_hdr = wcs.WCS(self.data_cube.header).celestial.to_header()
        hdulist.append(fits.ImageHDU(data=twodmap, header=wcs_hdr))
        if uncert_idx is not None:
            hdulist.append(fits.ImageHDU(data=twodmap_uncert, header=wcs_hdr))
        hdulist.writeto(outfile, clobber=clobber)

        print("2D map saved to {}".format(outfile))

    def make_emission_line_cube(self, clobber=False):
        """
        Subtract starlight continuum fits from data and save the cube.

        Where binning has been performed, the individual spectra are replaced
        by the binned spectrum. The emission line cube is saved with suffix
        `_emissionline.fits`.

        Parameters
        ----------
        clobber : bool, optional
            Whether to overwrite an existing emission line cube.
        """
        outfile = self.base_name+"_emissionline.fits"
        if not clobber and os.path.isfile(outfile):
                print("{} exists and clobber is false".format(outfile))
                return

        # make a dummy cube to hold the emission line spectra
        emission_line_cube = np.full(self.data_cube.shape, np.nan)
        for bn in self._get_bin_nums("nocontbad"):
            bin_res = self.results["bin"][bn]
            if bin_res["continuum"]["bad"] == 1:
                continue
            # subtract the synthetic spectrum from the data and account for
            # normalisation in starlight
            emission_line_spec = ((bin_res["continuum"]["sl_spec"][:, 1]
                                   - bin_res["continuum"]["sl_spec"][:, 2])
                                  * bin_res["continuum"]["fobs_norm"])
            for x, y in zip(*bin_res["spax"]):
                emission_line_cube[:, y, x] = emission_line_spec

        hdulist = fits.HDUList()
        hdulist.append(fits.ImageHDU(data=emission_line_cube,
                                     header=self.data_cube.header))
        hdulist.append(fits.ImageHDU(data=self.stddev_cube.data,
                                     header=self.stddev_cube.header))
        hdulist.writeto(outfile, clobber=clobber)

        print("emission line cube saved to {}".format(outfile))

    def make_continuum_cube(self):
        """
        Subtract the emission line model from each bin to produce a continuum
        cube.

        To be implemented.
        """
        pass
        # print("continuum cube saved to {}".format(outfile))

    def plot_continuum(self, bin_num):
        """
        Plot the spectrum and results of starlight continuum fitting for
        ``bin_num``
        """
        x_mean, y_mean = self.results["bin"][bin_num]["mean"]
        bin_res = self.results["bin"][bin_num]["continuum"]

        if bin_res["bad"]:
            print("bad continuum fitting for bin {}, skipping plotting"
                  .format(bin_num))
            return

        # Get the data we want to plot:
        lamb = bin_res["sl_spec"][:, 0]
        obs = np.ma.masked_array(bin_res["sl_spec"][:, 1],
                                 mask=bin_res["sl_spec"][:, 1] < 0)
        syn = bin_res["sl_spec"][:, 2]
        resid = obs - syn
        err = np.ma.masked_array(1/bin_res["sl_spec"][:, 3],
                                 mask=bin_res["sl_spec"][:, 3] < 0)
        masked = np.ma.masked_array(resid, mask=bin_res["sl_spec"][:, 3] == 0)
        slice_idxs = [(s.start, s.stop-1) for s in np.ma.clump_masked(masked)]
        clipped = np.ma.masked_array(resid, mask=bin_res["sl_spec"][:, 3] > -1)

        b = bin_res["bases"]
        Z = np.sort(np.unique(b[:, 5]))
        zages = []
        lightweights = []
        massweights = []
        for _Z in Z:
            _b = b[b[:, 5] == _Z]
            zages.append(np.log10(_b[:, 4]))
            lightweights.append(_b[:, 1])  # x_j
            massweights.append(_b[:, 3])  # M_cor #FIXME M_ini?????

        plt.close("all")
        # Plot the spectrum and continuum fit
        slfig = plt.figure()
        gs = gridspec.GridSpec(4, 2, height_ratios=[1, 1, 1, 1],
                               width_ratios=[3, 1])
        axol = slfig.add_subplot(gs[:3, 0])
        axol.plot(lamb, obs, c=OBSCOL, ls="-", lw=1, label="observed")
        axol.plot(lamb, syn, c=SYNTHCOL, ls="--", lw=1, label="starlight fit")
        axol.plot(lamb, err, c=ERRCOL, ls="-", lw=1, label="error")
        for idx in slice_idxs:
            axol.axvspan(self.lamb[idx[0]], self.lamb[idx[1]], color=MASKCOL,
                         alpha=0.3)
        axol.set_ylabel("F$_\\lambda$ [normalised]")
        axol.set_ylim(0, 2.2)
        plt.setp(axol.get_xticklabels(), visible=False)  # we share the x axis

        axresid = slfig.add_subplot(gs[3, 0], sharex=axol)
        axresid.plot(lamb, resid, c=OBSCOL, ls="-", lw=1, label="residual")
        axresid.plot(lamb, clipped, c=CLIPPEDCOL,
                     ls="-", lw=1, label="clipped")
        for idx in slice_idxs:
            axresid.axvspan(self.lamb[idx[0]], self.lamb[idx[1]],
                            color=MASKCOL, alpha=0.3)
        axresid.set_xlim(np.min(lamb), np.max(lamb))
        axresid.set_xlabel("Wavelength [{}]"
                           .format(self.lamb_unit.to_string("latex")))
        axresid.set_ylabel("Residual")
        axresid.set_yticks([-0.2, 0.0, 0.2, 0.4, 0.6])
        axresid.set_ylim(-0.3, 0.7)

        # Plot the contribution of the SSPs to the population
        axlight = slfig.add_subplot(gs[:2, 1])
        axmass = slfig.add_subplot(gs[2:, 1], sharex=axlight)
        bottoml, bottomm = np.zeros(len(zages[0])), np.zeros(len(zages[0]))
        for i, _Z in enumerate(Z):
            axlight.bar(zages[i], lightweights[i], color=ZCOL[i],
                        bottom=bottoml, label="Z$_\\star$ = "+str(_Z),
                        edgecolor="none", align="center", width=0.2)
            axmass.bar(zages[i], massweights[i], color=ZCOL[i], bottom=bottomm,
                       label="Z$_\\star$ = "+str(_Z), edgecolor="none",
                       align="center", width=0.2)
            bottoml += lightweights[i]
            bottomm += massweights[i]
        axlight.plot(zages[0], [77]*len(zages[0]), marker="|", markersize=12,
                     color=SYNTHCOL, ls="none")
        lgnd = axlight.legend(loc=9, bbox_to_anchor=(0.5, 1.23),
                              frameon=False, fontsize=10, ncol=2)
        axlight.set_ylabel("x$_\\textrm{{j}}$ [\%]")
        axlight.set_ylim(0, 80)
        plt.setp(axlight.get_xticklabels(), visible=False)  # share the x axis

        axmass.set_yscale("log")
        axmass.set_ylim(0.1, 105)
        yformatter = ticker.FuncFormatter(
            lambda y, pos: ('{{:.{:1d}f}}'.format(
                int(np.maximum(-np.log10(y), 0)))).format(y))
        axmass.yaxis.set_major_formatter(yformatter)
        axmass.set_xticks([6, 7, 8, 9, 10])
        axmass.set_xlim(5.5, 10.5)
        axmass.set_xlabel("$\\log_{10}$ Age$_\\star$ [years]")
        axmass.set_ylabel("Mcor$_\\textrm{{j}}$ [\%]")

        slfig.text(0.15, 1.0, "$\chi^{{2}}/\\textrm{{dof}} = {:.3f}$, "
                   "$\\bar{{x}} = {:.2f}$, $\\bar{{y}} = {:.2f}$"
                   .format(bin_res["chi2/Nl_eff"], x_mean, y_mean),
                   color="k", size=12)

        slfig.text(0.15, 0.96, "A$_\\textrm{{v}} = {:5.3f}$, "
                   "$\\sigma_\\star = {:6.2f}$ km s$^{{-1}}$, "
                   "$v_\\star = {:6.2f}$ km s$^{{-1}}$".format(
                       bin_res["AV_min"], bin_res["vd_min"],
                       bin_res["v0_min"]), color="k", size=12)
        slfig.tight_layout()
        slfig.subplots_adjust(hspace=0.1)
        slfig.savefig(self.base_name+"_sl_fit_{}.png".format(bin_num),
                      bbox_inches="tight", dpi=300, additional_artists=(lgnd,))
        print("plot saved to {}_sl_fit_{}.png".format(self.base_name, bin_num))

    def plot_yio(self, age1=5e8, age2=5e9):
        """
        Plot the contribution of the young, intermediate and old stellar
        populations.

        The plotted values are based on the light fraction results of the
        starlight continuum fitting. Nan or negative values are shown as white.
        The plot is saved with suffix `_yio.pdf`.

        Parameters
        ----------
        age1, age2: int, optional
            The dividing ages in years between the young (<= ``age1``),
            intermediate (between ``age1`` and ``age2``) and old (>= ``age2``)
            stellar components.
        """

        young = np.full(self.data_cube.shape[1:], np.nan)
        inter = np.full(self.data_cube.shape[1:], np.nan)
        old = np.full(self.data_cube.shape[1:], np.nan)
        age1_str = (str(int(age1/1e6)) + " Myr" if age1 < 1e9 else
                    str(int(age1/1e9)) + " Gyr")
        age2_str = (str(int(age2/1e6)) + " Myr" if age2 < 1e9 else
                    str(int(age2/1e9)) + " Gyr")
        for bn in self._get_bin_nums("nocontbad"):
            bin_res = self.results["bin"][bn]
            yngidx = bin_res["continuum"]["bases"][:, 4] <= age1
            intidx = ((age1 < bin_res["continuum"]["bases"][:, 4])
                      & (bin_res["continuum"]["bases"][:, 4] < age2))
            oldidx = bin_res["continuum"]["bases"][:, 4] >= age2
            ynglightfrac = np.sum(bin_res["continuum"]["bases"][yngidx, 1])
            intlightfrac = np.sum(bin_res["continuum"]["bases"][intidx, 1])
            oldlightfrac = np.sum(bin_res["continuum"]["bases"][oldidx, 1])
            young[bin_res["spax"][::-1]] = ynglightfrac
            inter[bin_res["spax"][::-1]] = intlightfrac
            old[bin_res["spax"][::-1]] = oldlightfrac

        plt.close("all")
        axyoung = plt.subplot(311, adjustable="box-forced")
        plt.imshow(young, origin="lower", interpolation="none", cmap="Blues",
                   vmin=0, vmax=np.nanmax(young))
        axyoung.autoscale(False)
        plt.plot(self.nucleus[0], self.nucleus[1], "kx", markersize=10)
        plt.title("Age$_\\star <$ {}".format(age1_str))
        plt.colorbar(label="x$_\\textrm{j}$ [\%]")

        axinter = plt.subplot(312, sharex=axyoung, sharey=axyoung,
                              adjustable="box-forced")
        plt.imshow(inter, origin="lower", interpolation="none", cmap="Greens",
                   vmin=0, vmax=np.nanmax(inter))
        axinter.autoscale(False)
        plt.plot(self.nucleus[0], self.nucleus[1], "kx", markersize=10)
        plt.title("{} $<$ Age$_\\star <$ {}".format(age1_str, age2_str))
        plt.colorbar(label="x$_\\textrm{j}$ [\%]")

        axold = plt.subplot(313, sharex=axyoung, sharey=axyoung,
                            adjustable="box-forced")
        plt.imshow(old, origin="lower", interpolation="none", cmap="Reds",
                   vmin=0, vmax=np.nanmax(old))
        axold.autoscale(False)
        plt.plot(self.nucleus[0], self.nucleus[1], "kx", markersize=10)
        plt.title("Age$_\\star >$ {}".format(age2_str))
        plt.colorbar(label="x$_\\textrm{j}$ [\%]")

        fig = plt.gcf()
        fig.set_size_inches(5, 15)

        plt.savefig(self.base_name+"_yio.pdf", bbox_inches="tight")
        print("plot saved to {}".format(self.base_name+"_yio.pdf"))

    def plot_kinematics(self, norm_v0=0):
        """
        Plot the kinematics of the host.

        The plotted values are based on the velocity offset and dispersion
        results of the starlight continuum fitting. Nan or negative values are
        shown as white. The plot is saved with suffix `_kinematics.pdf`.

        Parameters
        ----------
        norm_v0: float or str, optional
            The zero-point of the velocity distribution (i.e. the plotted
            values will be v0 - norm_v0). The special case "nucleus" will set
            the zero-point as the value of the nucelus bin.
        """

        v0 = np.full(self.data_cube.shape[1:], np.nan)
        vd = np.full(self.data_cube.shape[1:], np.nan)
        if norm_v0 == "nucleus":
            bn = self.get_loc_bin(self.nucleus)
            if bn is None:
                print("nucleus is not in a bin")
                return
            norm_v0 = self.results["bin"][bn]["continuum"]["v0_min"]
        elif not isinstance(norm_v0, (float, int)):
            print("norm_v0 must be 'nucleus' or a float/int")
            return
        for bn in self._get_bin_nums("nocontbad"):
            bin_res = self.results["bin"][bn]
            v0[bin_res["spax"][::-1]] = (bin_res["continuum"]["v0_min"]
                                         - norm_v0)
            vd[bin_res["spax"][::-1]] = bin_res["continuum"]["vd_min"]

        plt.close("all")
        v0min, v0max = np.nanmin(v0), np.nanmax(v0)
        v0cmap = shiftedColorMap(cm.coolwarm,
                                 midpoint=(1. - v0max/(v0max + abs(v0min))))
        ax1 = plt.subplot(121, adjustable="box-forced")
        plt.imshow(v0, origin="lower", interpolation="none", cmap=v0cmap)
        plt.colorbar(label="$v_\\star$ [km~s$^{-1}$]",
                     orientation="horizontal").ax.tick_params(labelsize=10)
        ax1.autoscale(False)
        plt.plot(self.nucleus[0], self.nucleus[1], "kx", markersize=10)

        ax2 = plt.subplot(122, sharex=ax1, sharey=ax1, adjustable="box-forced")
        plt.imshow(vd, origin="lower", interpolation="none", cmap="afmhot_r")
        plt.colorbar(label="$\\sigma_\\star$ [km~s$^{-1}$]",
                     orientation="horizontal").ax.tick_params(labelsize=10)
        ax2.autoscale(False)
        plt.plot(self.nucleus[0], self.nucleus[1], "kx", markersize=10)

        plt.savefig(self.base_name+"_kinematics.pdf", bbox_inches="tight")
        print("plot saved to {}".format(self.base_name+"_kinematics.pdf"))

    def plot_emission(self, bin_num):
        """
        Plot the residual emission line spectrum and results of emission
        line fitting for ``bin_num``

        """
        bin_res = self.results["bin"][bin_num]

        lamb = bin_res["continuum"]["sl_spec"][:, 0]
        mask = bin_res["spec"][:, 3] == 2

        emline_obs = np.ma.masked_array(
            (bin_res["continuum"]["sl_spec"][:, 1] -
             bin_res["continuum"]["sl_spec"][:, 2]) *
            bin_res["continuum"]["fobs_norm"], mask=mask)
        emline_uncert = np.ma.masked_array(bin_res["spec"][:, 2], mask=mask)
        emline_model = _get_emline_model(self._emission_lines,
                                         bin_res["emission"])
        # The additional residual function that was subtracted from the
        # model. This will be 0.0 if not used so not affect the emline spec
        resid_fn = bin_res["emission"]["resid_fn"]
        emline_obs -= resid_fn

        plt.close("all")
        # Determine how many zoomed-in windows on the full wavelength
        # range we need
        lambm = np.ma.masked_array(lamb, mask=np.ones(len(lamb)))
        for submodel in emline_model:
            low_idx = np.abs(lamb - (submodel.mean-25)).argmin()
            upp_idx = np.abs(lamb - (submodel.mean+25)).argmin()
            lambm.mask[low_idx:upp_idx] = 0
        slices = np.ma.clump_unmasked(lambm)
        nax = len(slices)
        elfig, axes = plt.subplots(1, nax, sharey=True, figsize=(13, 5))
        if nax == 1:
            axes = [axes]
        for slc, ax in zip(slices, axes):
            ax.plot(lamb[slc], emline_obs[slc], c=OBSCOL, lw=2)
            ax.fill_between(lamb[slc], emline_obs[slc]-emline_uncert[slc],
                            emline_obs[slc]+emline_uncert[slc], color=ERRCOL)
            x = np.arange(lamb[slc.start], lamb[slc.stop], 0.25)
            ax.plot(x, emline_model(x), c=MASKCOL, lw=2, ls="--")
            for submodel in emline_model:
                if not (lamb[slc.start] < submodel.mean < lamb[slc.stop]):
                    continue
                ax.axvline(submodel.mean.value, ls="--", c="grey", lw=0.5)
                label = " ".join(submodel.name.split("_")[:2])
                ax.text(submodel.mean.value-2, submodel.amplitude.value*0.7,
                        label, rotation=90, ha="right", va="bottom")
            ax.set_xlim(lamb[slc.start], lamb[slc.stop])
            ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
            ax.tick_params(axis="x", labelsize=10)

        axes[0].spines["right"].set_visible(False)
        axes[-1].spines["left"].set_visible(False)
        axes[0].yaxis.tick_left()
        axes[-1].yaxis.tick_right()
        for ax in axes[1:-1]:
            ax.yaxis.set_ticks_position("none")
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)

        # Add common x and y labels
        bigax = elfig.add_subplot(111, frameon=False)
        bigax.set_facecolor("none")
        bigax.tick_params(labelcolor="none", top="off", bottom="off",
                          left="off", right="off")
        bigax.set_ylabel("Flux [{}]".format(self.flux_unit))
        bigax.set_xlabel("Wavelength [{}]"
                         .format(self.lamb_unit.to_string("latex")))

        elfig.suptitle("$\\chi^2/\\textrm{{dof}} = {:.3f}$, " "$\\bar{{x}} = "
                       "{:.2f}$, $\\bar{{y}} = {:.2f}$"
                       .format(bin_res["emission"]["chi2dof"],
                               bin_res["mean"][0], bin_res["mean"][1]),
                       color="k", size=12)
        elfig.subplots_adjust(wspace=0.1)
        elfig.tight_layout()
        elfig.savefig(self.base_name+"_el_fit_{}.png".format(bin_num),
                      bbox_inches="tight", dpi=300)
        print("plot saved to {}_el_fit_{}.png".format(self.base_name, bin_num))

    def plot_metallicity(self, indicator="D16", cumweight=False):
        """
        Plot the metallicity of the host.

        A map of the host is shown with bins colour coded by their metallicity
        in the chosen indicator. Nan or negative values are shown as white. The
        plot is saved with suffix `_metallicity_`\ ``indicator``\ `.pdf`.

        Parameters
        ----------
        indicator: str
            The metallicity indicator to plot, current options are "PP04_N2",
            "PP04_O3N2", "M13_N2", "M13_O3N2", "D16".
        cumweight: bool, optional
            If ``True``, the contribution of each bin to the cumulative
            histogram of metallicities will be weighted by its flux in
            "Halpha_6563".
        """
        valid_indicators = ["PP04_N2", "PP04_O3N2",
                            "M13_N2", "M13_O3N2",
                            "D16"]
        if indicator not in valid_indicators:
            raise AttributeError("`indicator` must be one of {}"
                                 .format(valid_indicators))

        bin_nums, n_custom = self._get_bin_nums("nobad", custom=True)
        # Initialise an empty map to populate with Z values
        Zmap = np.full(self.data_cube.shape[1:], np.nan)
        # Store the metallicity values and uncertainties
        Zvals = np.full((len(bin_nums), 2), np.nan)
        # and Halpha flux for cumulative weighting, if needed
        Wvals = np.full((len(bin_nums)), np.nan)
        # Hold the mean, min and max distance of the bins from the nucleus
        dists = np.full((len(bin_nums), 3), np.nan)

        for i, bn in enumerate(bin_nums):
            bin_res = self.results["bin"][bn]
            Z = bin_res["emission"]["metallicity"][indicator]
            Zmap[bin_res["spax"][::-1]] = Z[0]
            Zvals[i] = Z
            dists[i] = (bin_res["dist_mean"], bin_res["dist_min"],
                        bin_res["dist_max"])
            if cumweight:
                Wvals[i] = (bin_res["emission"]["lines"]["Halpha_6563"]
                            ["flux"][0])
        # Remove any bins without a metallicity determined
        Zvalsnn = Zvals[~np.isnan(Zvals[:, 0]), 0]
        if cumweight:
            weights = Wvals[~np.isnan(Zvals[:, 0])]
        else:
            weights = None

        plt.close("all")
        zfig = plt.figure()
        gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1],
                               width_ratios=[2, 1])
        # A map of the bins and their Z values
        axmap = zfig.add_subplot(gs[:, 0])
        m = axmap.imshow(Zmap, origin="lower", interpolation="none")
        c = plt.colorbar(mappable=m, ax=axmap, orientation="horizontal",
                         label="$Z$ [$12 + \log_{10}(\\textrm{O/H})$]")
        c.ax.tick_params(labelsize=16)
        axmap.autoscale(False)
        axmap.plot(self.nucleus[0], self.nucleus[1], "kx", markersize=10)

        # A cumulative histogram of the Z values
        axcum = zfig.add_subplot(gs[0, 1])
        bin_edges = np.append(np.sort(Zvalsnn), np.max(Zvalsnn)+1e-9)
        n, b, p = axcum.hist(Zvalsnn, weights=weights, cumulative=True,
                             normed=True, histtype="step", linewidth=3,
                             bins=bin_edges, color=c.cmap(0.3))
        #    and show the custom bins highlighted
        for i in range(n_custom):
            axcum.axvline(Zvals[i, 0], color=c.cmap(0.9), lw=2)
        p[0].set_xy(p[0].get_xy()[:-1])
        axcum.set_xlim(min(Zvalsnn), np.max(Zvalsnn))
        axcum.tick_params(axis="x", labelsize=10)
        axcum.set_ylim(0, 1)
        axcum.set_xlabel("$Z$ [$12 + \log_{10}(\\textrm{O/H})$]")
        axcum.set_ylabel("Cumulative Fraction")

        # A plot of the Z value vs radial distance from nucleus
        axrad = zfig.add_subplot(gs[1, 1])
        axrad.errorbar(dists[:, 0], Zvals[:, 0],
                       xerr=[(dists[:, 0]-dists[:, 1]),
                       (dists[:, 2]-dists[:, 0])], yerr=Zvals[:, 1],
                       color=c.cmap(0.3), ls="none", marker="o", mew=0.3,
                       ms=4, capsize=0, ecolor=c.cmap(0.3))
        #    and show the custom bins highlighted
        axrad.errorbar(dists[:n_custom, 0], Zvals[:n_custom, 0],
                       xerr=[(dists[:n_custom, 0]-dists[:n_custom, 1]),
                       (dists[:n_custom, 2]-dists[:n_custom, 0])],
                       yerr=Zvals[:n_custom, 1],
                       color=c.cmap(0.9), ls="none", marker="*", mew=0.3,
                       ms=9, capsize=0, ecolor=c.cmap(0.9))
        axrad.set_xlabel("Distance from nucleus [pixels]")
        axrad.set_ylabel("$Z$ [$12 + \log_{10}(\\textrm{O/H})$]")
        zfig.suptitle("Metallicity indicator {}"
                      .format(indicator.replace("_", " ")))
        gs.tight_layout(zfig, rect=[0, 0.03, 1, 0.97])
        zfig.savefig(self.base_name+"_metallicity_{}.pdf".format(indicator),
                     bbox_inches="tight")
        print("plot saved to {}_metallicity_{}.pdf".format(self.base_name,
                                                           indicator))

    def plot_bpt(self, snrlimit=3):
        """
        Plot the BPT diagram for all bins

        The [NII]/H\ :math:`\\alpha` vs [OIII]/H\ :math:`\\beta` [BPT]_ diagram
        is plotted for each bin, along with a 2D map showing the classification
        of each bin. Classification dividing relations are taken from [K13]_
        for the AGN-HII division, and [K01]_ for the theretical star-formation
        driven limit.

        Parameters
        ----------
        snrlimit: float, optional
            Only bins where all lines were detected above this value will be
            plotted.

        References
        ----------
        .. [BPT] Baldwin, J. A., Phillips, M. M. & Terlevich, R.,
                 "Classification parameters for the emission-line spectra of
                 extragalactic objects", 1981, PASP, 93, 53
        .. [K01] Kewley, L. J. et al. "Theoretical Modeling of Starburst
                 Galaxies", ApJ, 556, 121
        .. [K13] Kewley, L. J. et al. "The Cosmic BPT Diagram: Confronting
                 Theory with Observations", 2013, ApJL, 774, 10
        """
        # Redshift dependant dividing relation between Hii regions and AGN
        # from Kewley et al. 2013 eq. 1
        # x = [NII]/Halpha flux ratio
        OIIIHb_k13 = lambda x, z=self._z: (0.61/(np.log10(x)-0.02-0.1833*z)
                                           + 1.2 + 0.03 * z)
        # and from Kewley et al. 2001 eq. 5 the maximal starformation limit
        OIIIHb_k01 = lambda x: (0.61/(np.log10(x)-0.47) + 1.19)

        bin_nums, n_custom = self._get_bin_nums("nobad", custom=True)
        # Initialise an empty map to populate with values denoting:
        # 0 = Hii, 1 = Hii-> maximal starburst, 2 = AGN
        valmap = np.full(self.data_cube.shape[1:], np.nan)
        # Hold the line ratios for each bin:
        # ([NII]/Ha, uncert, [OIII]/Hb, uncert, val) per row
        lineratios = np.full((len(bin_nums), 5), np.nan)

        lines = ("[NII]_6583", "Halpha_6563", "[OIII]_5007", "Hbeta_4861")
        for i, bn in enumerate(bin_nums):
            bin_res = self.results["bin"][bn]
            snr = np.array([bin_res["emission"]["lines"][line]["snr"]
                            for line in lines])
            # Do not plot for any bins below snrlimit
            if not all(snr > snrlimit):
                continue
            f = [bin_res["emission"]["lines"][line]["flux"][0] for line in
                 lines]
            fu = [bin_res["emission"]["lines"][line]["flux"][1] for line in
                  lines]
            NIIHa = f[0]/f[1]
            NIIHa_uncert = NIIHa * ((fu[0]/f[0])**2 + (fu[1]/f[1])**2)**0.5
            OIIIHb = f[2]/f[3]
            OIIIHb_uncert = OIIIHb * ((fu[2]/f[2])**2 + (fu[3]/f[3])**2)**0.5
            val = 0  # The value to show as on map
            if math.log10(OIIIHb) > OIIIHb_k13(NIIHa):
                val += 1  # above HII/AGN divider
            if math.log10(OIIIHb) > OIIIHb_k01(NIIHa):
                val += 1  # above maximal starburst
            valmap[bin_res["spax"][::-1]] = val
            lineratios[i] = [math.log10(NIIHa), 0.434*(NIIHa_uncert/NIIHa),
                             math.log10(OIIIHb), 0.434*(OIIIHb_uncert/OIIIHb),
                             val]

        plt.close("all")
        vfig = plt.figure()
        gs = gridspec.GridSpec(1, 2)
        # A map of the bins and their values
        axmap = vfig.add_subplot(gs[0, 0])
        c = cm.viridis
        cmap = colors.ListedColormap([c(0.1), c(0.5), c(0.9)])
        norm = colors.BoundaryNorm([0, 1, 2, 3], cmap.N)
        m = axmap.imshow(valmap, origin="lower", interpolation="none",
                         cmap=cmap, norm=norm)
        cbar = plt.colorbar(mappable=m, ax=axmap, orientation="horizontal")
        cbar.ax.get_xaxis().set_ticks([])
        for j, lab in enumerate(["HII", "COMP", "AGN"]):
            cbar.ax.text((2*j + 1) / 6., -1.5, lab, ha="center")
        cbar.ax.tick_params(labelsize=16)
        axmap.autoscale(False)
        axmap.plot(self.nucleus[0], self.nucleus[1], "kx", markersize=10)
        # The BPT scatter plot
        # get the rows for each type so we can colour differently on plot
        hii = lineratios[:, 4] == 0
        sb = lineratios[:, 4] == 1
        agn = lineratios[:, 4] == 2
        axrat = vfig.add_subplot(gs[0, 1])
        x0 = np.linspace(-2, 0.02+0.1832*self._z, 100)
        x1 = np.linspace(-2, 0.46, 100)
        axrat.plot(x0, OIIIHb_k13(10**x0), "k-")
        axrat.plot(x1, OIIIHb_k01(10**x1), "k--", label="Kewley+01")
        for r, cval in zip((hii, sb, agn), (0.1, 0.5, 0.9)):
            axrat.errorbar(lineratios[r, 0], lineratios[r, 2],
                           xerr=lineratios[r, 1], yerr=lineratios[r, 3],
                           color=c(cval), ls="none", marker="o", mew=0.3,
                           ms=6, capsize=0, ecolor=c(cval))
        axrat.legend(loc=0, frameon=False)
        axrat.set_xlabel("$\\log_{10}([\\textrm{NII}]/\\textrm{H}\\alpha)$")
        axrat.set_ylabel("$\\log_{10}([\\textrm{OIII}]/\\textrm{H}\\beta)$")
        axrat.set_ylim(min(-1, np.min(lineratios[:, 1])), 1.5)
        axrat.text(-1.6, -0.2, "\\bf{HII}", size=20, color="darkgrey")
        axrat.text(-0.15, 0.85, "\\bf{AGN}", size=20, color="darkgrey")
        vfig.tight_layout()
        vfig.savefig(self.base_name+"_bpt.pdf", bbox_inches="tight")
        print("plot saved to {}_bpt.pdf".format(self.base_name))

    def plot_line_map(self, line, snrlimit=3):
        """
        Plot maps of emission line parameters.

        Plots 2D maps and radial distributions of the EW, flux, velocity offset
        and FWHM of the chosen line.

        Parameters
        ----------
        line: str
            The emission line to plot. This will take the form of,
            e.g. "Halpha_6563", "[NII]_6583". i.e the line name and the rounded
            wavelength of the line as per the data in the emission lines json
            data given to :class:`~ifuanal.IFUCube`. Set as "?" to see a list
            available.
        snrlimit: float, optional
            Only bins where the line was detected above this value will be
            plotted.
        """
        for bn in self._get_bin_nums("nobad"):
            try:
                lines = self.results["bin"][bn]["emission"]["lines"].keys()
            except KeyError:
                pass
            else:
                if line not in lines or line == "?":
                    print("available emission lines:")
                    print(", ".join(lines))
                    return
                else:
                    break

        bin_nums, n_custom = self._get_bin_nums("nobad", custom=True)
        # Initialise an empty 4-layer map to populate
        val_maps_shape = (4, self.data_cube.shape[1], self.data_cube.shape[2])
        val_maps = np.full(val_maps_shape, np.nan)
        # Also hold the data and uncerts in a list for the radial plots
        val_list = np.full((len(bin_nums), 4), np.nan)
        val_uncert_list = np.full((len(bin_nums), 4), np.nan)
        # Hold the mean, min and max distance of the bins from the nucleus
        dists = np.full((len(bin_nums), 3), np.nan)

        for i, bn in enumerate(bin_nums):
            bin_res = self.results["bin"][bn]
            # Only consider those bins with SNR above our limit
            if bin_res["emission"]["lines"][line]["snr"] < snrlimit:
                continue
            for j, prop in enumerate(("ew", "flux", "offset", "fwhm")):
                val, val_uncert = bin_res["emission"]["lines"][line][prop]
                val_maps[j][bin_res["spax"][::-1]] = val
                val_list[i, j] = val
                val_uncert_list[i, j] = val_uncert
                dists[i] = (bin_res["dist_mean"], bin_res["dist_min"],
                            bin_res["dist_max"])
        plt.close("all")
        lfig = plt.figure(figsize=(16, 8))
        gs = gridspec.GridSpec(2, 4, height_ratios=[2, 1])
        propdict = {0: {"prop": "ew",
                        "name": "Equivalent Width",
                        "unit": self.lamb_unit.to_string("latex")},
                    1: {"prop": "flux",
                        "name": "Flux",
                        "unit": self.flux_unit},
                    2: {"prop": "offset",
                        "name": "Velocity offset",
                        "unit": "km s$^{{-1}}$"},
                    3: {"prop": "fwhm",
                        "name": "FWHM",
                        "unit": "km s$^{{-1}}$"}
                    }
        for i, v in propdict.items():
            # Map of bins and their values for each prop
            axmap = lfig.add_subplot(gs[0, i])
            m = axmap.imshow(val_maps[i], origin="lower", interpolation="none")
            c = plt.colorbar(mappable=m, ax=axmap, orientation="horizontal",
                             label="{} [{}]".format(v["name"], v["unit"]))
            c.ax.tick_params(labelsize=12)
            axmap.autoscale(False)
            axmap.plot(self.nucleus[0], self.nucleus[1], "kx", markersize=10)
            axmap.plot(self.nucleus[0], self.nucleus[1], "kx", markersize=10)

            # A radial plot below the map
            axrad = lfig.add_subplot(gs[1, i])
            axrad.errorbar(dists[:, 0], val_list[:, i],
                           xerr=[(dists[:, 0]-dists[:, 1]),
                                 (dists[:, 2]-dists[:, 0])],
                           yerr=val_uncert_list[:, i], color=c.cmap(0.3),
                           ls="none", marker="o", mew=0.3, ms=4, capsize=0,
                           ecolor=c.cmap(0.3))
            # and show the custom bins highlighted
            axrad.errorbar(dists[:n_custom, 0], val_list[:n_custom, i],
                           xerr=[(dists[:n_custom, 0]-dists[:n_custom, 1]),
                                 (dists[:n_custom, 2]-dists[:n_custom, 0])],
                           yerr=val_uncert_list[:n_custom, i],
                           color=c.cmap(0.9), ls="none", marker="*", mew=0.3,
                           ms=9, capsize=0, ecolor=c.cmap(0.9))
            axrad.set_xlabel("Distance from nucleus [pixels]")
            axrad.set_ylabel("{} [{}]".format(v["name"], v["unit"]))
        lfig.suptitle("Emission line {}".format(line.replace("_", " ")))
        gs.tight_layout(lfig, rect=[0, 0.03, 1, 0.97])
        lfig.savefig(self.base_name+"_{}.pdf".format(line),
                     bbox_inches="tight")
        print("plot saved to {}_{}.pdf".format(self.base_name, line))

    def plot_extinction(self, ebv=False):
        """
        Plot the extinction map of the host.

        The plotted values are based on the Av_min parameter found by the
        starlight continuum fitting for the stellar extinction and via the
        Balmer decrement for the ionised gas component. Nan values are
        shown as white. The plot is saved with suffix `_extinction.pdf`.

        Parameters
        ----------
        ebv: bool, optional
            If True will plot the values as E(B-V) values, using ``RV``,
            otherwise AV is shown
        """

        label = ("$E(B-V)_{}$" if ebv
                 else "A$_{{V, {}}}$")
        c = 1. if ebv else self.RV

        stellar = np.full(self.data_cube.shape[1:], np.nan)
        gas = np.full(self.data_cube.shape[1:], np.nan)
        for bn in self._get_bin_nums("nocontbad"):
            bin_res = self.results["bin"][bn]
            stellar[bin_res["spax"][::-1]] = bin_res["continuum"]["ebv_star"]*c
            gas[bin_res["spax"][::-1]] = bin_res["emission"]["ebv_gas"]*c

        plt.close("all")
        smin, smax = np.nanmin(stellar), np.nanmax(stellar)
        if smin < 0:
            stellarcmap = shiftedColorMap(cm.coolwarm,
                                          midpoint=(1.
                                                    - smax/(smax + abs(smin))))
        else:
            stellarcmap = cm.Reds
        ax1 = plt.subplot(121, adjustable="box-forced")
        plt.imshow(stellar, origin="lower", interpolation="none",
                   cmap=stellarcmap)
        plt.colorbar(label=label.format("\\star"),
                     orientation="horizontal").ax.tick_params(labelsize=10)
        ax1.autoscale(False)
        plt.plot(self.nucleus[0], self.nucleus[1], "kx", markersize=10)

        gmin, gmax = np.nanmin(gas), np.nanmax(gas)
        if gmin < 0:
            gascmap = shiftedColorMap(cm.coolwarm,
                                      midpoint=(1. - gmax/(gmax + abs(gmin))))
        else:
            gascmap = cm.Reds
        ax2 = plt.subplot(122, sharex=ax1, sharey=ax1, adjustable="box-forced")
        plt.imshow(gas, origin="lower", interpolation="none", cmap=gascmap)
        plt.colorbar(label=label.format("\\textrm{gas}"),
                     orientation="horizontal").ax.tick_params(labelsize=10)
        ax2.autoscale(False)
        plt.plot(self.nucleus[0], self.nucleus[1], "kx", markersize=10)

        plt.savefig(self.base_name+"_extinction.pdf", bbox_inches="tight")
        print("plot saved to {}".format(self.base_name+"_extinction.pdf"))

    def plot_worst_fits(self, N=5):
        """
        Find the N bins with the worst chi2 from the starlight and emission
        line fits and plot both fits for these bins.

        """
        c = self._get_bin_nums("nocontbad")
        worst_sl_bins = sorted(c, key=(lambda key:
                                       self.results["bin"][key]["continuum"]
                                       .get("chi2/Nl_eff", -1)))[-N:]
        print("worst continuum fitted bins: {}".format(worst_sl_bins))
        e = self._get_bin_nums("nobad")
        worst_el_bins = sorted(e, key=(lambda key:
                                       self.results["bin"][key]["emission"]
                                       ["chi2dof"]))[-N:]
        print("worst emission line fitted bins: {}".format(worst_el_bins))
        worst_bins = set(worst_sl_bins + worst_el_bins)
        for wb in worst_bins:
            self.plot_continuum(wb)
            self.plot_emission(wb)

    @classmethod
    def load_pkl(self, pkl_file):
        """
        Load a previous instance of IFUCube from its pickle file.

        Parameters
        ----------
        pkl_file : str
            The filepath of the pickled instance
        """
        try:
            cls = IFUCube.__new__(IFUCube)
            with open(pkl_file, "rb") as pkl_data:
                cls.__dict__ = pickle.load(pkl_data)
        except IOError:
            raise IOError("Couldn't create instance from pickle file"
                          " {}".format(pkl_file))
        # Read the cube HDUs from the saved FITS file
        cube_hdu = fits.open(pkl_file+".fits")
        cls.__dict__["prim_cube"] = cube_hdu[0]
        cls.__dict__["data_cube"] = cube_hdu[1]
        cls.__dict__["stddev_cube"] = cube_hdu[2]
        # Let's update n_cpu incase we're loading on a different machine
        cls.n_cpu = int(min(mp.cpu_count()-1, mp.cpu_count()*0.9))
        # Update the base_name in case we're in a different location
        cls.base_name = os.path.splitext(os.path.abspath(pkl_file))[0]
        print("loaded pkl file {}".format(pkl_file))
        return cls

    def save_pkl(self, pkl_file=None, clobber=False):
        """
        Load a previous instance of IFUCube from its pickle file.

        If the image was loaded from a compressed format (e.g. .fz), it will
        be stored in the pickle file as uncompressed as compressed files
        can't be pickled -- this will make the pickle file large, esp. for
        moasiacs!

        Parameters
        ----------
        pkl_file : None or str, optional
            The filepath to pickle the instance to. (defaults to ``None``, i.e.
            will add suffix ".pkl" to original datacube name)
        clobber : bool, optional
            Whether to overwrite an existing pkl and FITS file
        """
        if pkl_file is None:
            pkl_file = self.base_name+".pkl"

        try:
            open(pkl_file)
            open(pkl_file+".fits")
        except IOError:
            pass
        else:
            if not clobber:
                print("{0} and/or {0}.fits exist.".format(pkl_file))
                print("Use clobber=True to overwrite")
                return

        # Write the cube HDUs to a fits file as they may be large!
        tempcube = tempfile.mkstemp(prefix="ifuanal_", suffix=".pkl.fits",
                                    dir=os.path.dirname(self.base_name))[1]
        print("writing cube to temporary file {}".format(tempcube))
        cube_hdu = fits.HDUList([self.prim_cube,
                                 self.data_cube,
                                 self.stddev_cube])
        cube_hdu.writeto(tempcube)
        print("moving to {}".format(pkl_file+".fits"))
        shutil.move(tempcube, pkl_file+".fits")

        # Remove these entries from the dictionary we want to pickle
        pkl_dict = {}
        for key, val in self.__dict__.items():
            if key not in ["prim_cube", "data_cube", "stddev_cube"]:
                pkl_dict[key] = val

        # Write to a temporary file for safety - if there's pickling
        # errors, we don't want to overwrite our previous pickle.
        temp = tempfile.mkstemp(prefix="ifuanal_", suffix=".pkl",
                                dir=os.path.dirname(self.base_name))[1]
        print("writing instance to temporary pickle file {}".format(temp))
        with open(temp, "wb") as output:
            pickle.dump(pkl_dict, output)
        print("moving to {}".format(pkl_file))
        shutil.move(temp, pkl_file)


class MUSECube(IFUCube):
    """
    A Child class of :class:`~ifuanal.IFUCube` tailored for MUSE data cubes.

    Handles the conversion of the ERR extension (which is natively the
    variance) to std dev. Adds two headers used by :class:`~ifuanal.IFUCube`:
    IFU_EBV, IFU_Z.

    Parameters
    ----------
    muse_cube : str
        The filepath of the MUSE cube to analyse
    redshift : float
        The redshift of the observed host
    ebv : str or float, optional
        The value of Galactic extinction towards the hosts. The special
        case "IRSA" will contact the NASA IRSA service and automatically
        grab the value based on [SF11]_.
    RV : float, optional
        The RV value to use for the CCM extinction law.
    sl_dir : None or str, optional
        The directory containing starlight files and bases. The default
        ``None`` will use the `starlight/` subdirectory.
    el_json : None or str, optional
        Either a dictonary of emission lines to fit or the filepath to an
        emission lines json file. If ``None`` will use the default located in
        `data/emission_lines.json`. Follow the format of this default to create
        custom line lists.

    References
    ----------
    .. [SF11] Schlafly, E. & Finkbeiner, D., "Measuring Reddening with Sloan\
       Digital Sky Survey Stellar Spectra and Recalibrating SFD", ApJ, 2011

    """
    def __init__(self, muse_cube, redshift, ebv="IRSA", RV=3.1, sl_dir=None,
                 el_json=None):
        cube_hdu = fits.open(muse_cube)
        base_name = os.path.splitext(muse_cube)[0]
        if base_name.endswith(".fits"):
            # in case original file was, e.g., *.fits.fz
            base_name = base_name[:-5]
        # The MUSE STAT cube extension is the variance of the data, we want the
        # standard deviation
        cube_hdu[2].data = np.sqrt(cube_hdu[2].data)

        if ebv == "IRSA":
            ebv = get_IRSA_ebv(cube_hdu[1].header)
        cube_hdu[0].header["IFU_EBV"] = float(ebv)

        # Add in the redshift of the target since MUSE doesn't have this in the
        # header
        cube_hdu[0].header["IFU_Z"] = float(redshift)

        super(MUSECube, self).__init__(cube_hdu, base_name, RV,
                                       sl_dir, el_json)


def _get_emline_model(emlines, res=None):
    """
    Dynamically create the model to describe the emission lines.

    If res is given, the parameters of the model are filled from
    that dict.
    """
    model = None
    for elem, wl in emlines.items():
        for n, _wl in enumerate(map(float, wl)):
            submodel_name = "{}_{:.0f}_{}".format(elem, _wl, n)
            if res is not None:
                line = "{}_{:.0f}".format(elem, _wl)
                params = res["lines"][line]["fit_params"]
            else:
                params = (1., 0., 1.)
            if model is None:
                model = models.Gaussian1D(*params).rename(submodel_name)
            else:
                model += models.Gaussian1D(*params).rename(submodel_name)

    return model


def get_IRSA_ebv(header):
    """
    Retrieve the S&F 2011 Galactic E(B-V) for the reference pixel
    coordinates in a FITS header.
    """

    from astroquery.irsa_dust import IrsaDust
    coo = SkyCoord(header["CRVAL1"], header["CRVAL2"],
                   unit=(header["CUNIT1"], header["CUNIT2"]))
    tbl = IrsaDust.get_query_table(coo, section="ebv")
    ebv = float(tbl["ext SandF mean"])
    print("for location '{:}' found E(B-V) of {:5.3f}".format(coo.to_string(),
                                                              ebv))
    return ebv


def get_Alamb(lamb, ebv, RV=3.1, lamb_unit=u.Unit("angstrom")):
    """
    Return the wavelength dependant extinction for the input lambda array
    using the polynomial of [CCM89]_

    References
    ----------
    .. [CCM89] Cardelli, Clayton \& Mathis, 1989, ApJ, 345, 245, "The \
       relationship between infrared, optical, and ultraviolet extinction"
    """
    lamb = np.atleast_1d(lamb)
    # x in CCM is in 1/microns
    x = 1/(lamb * lamb_unit.to("micron"))

    a, b = np.zeros(x.shape, x.dtype), np.ndarray(x.shape, x.dtype)

    if any((x < 0.3) | (8 < x)):
            raise ValueError("Wavelengths extend beyond CCM extinction curve")

    ir = (0.3 <= x) & (x <= 1.1)
    opt = (1.1 <= x) & (x <= 3.3)
    nuv1 = (3.3 <= x) & (x <= 5.9)
    nuv2 = (5.9 <= x) & (x <= 8)

    # IR
    a[ir] = 0.574 * x[ir]**1.61
    b[ir] = -0.527 * x[ir]**1.61

    # opt/nIR
    a[opt] = np.polyval((0.32999, -0.7753, 0.01979, 0.72085, -0.02427,
                         -0.50447, 0.17699, 1), x[opt]-1.82)
    b[opt] = np.polyval((-2.09002, 5.3026, -0.62251, -5.38434, 1.07233,
                         2.28305, 1.41338, 0), x[opt]-1.82)

    # nUV1
    a[nuv1] = 1.752 - 0.316*x[nuv1] - 0.104 / ((x[nuv1] - 4.67)**2 + 0.341)
    b[nuv1] = -3.09 + 1.825*x[nuv1] + 1.206 / ((x[nuv1] - 4.62)**2 + 0.263)

    # nUV2
    y = x[nuv2] - 5.9
    Fa = -0.04473 * y**2 - 0.009779 * y**3
    Fb = -0.2130 * y**2 - 0.1207 * y**3
    a[nuv2] = (1.752 - 0.316 * x[nuv2] - 0.104
               / ((x[nuv2] - 4.67)**2 + 0.341) + Fa)
    b[nuv2] = (-3.09 + 1.825 * x[nuv2] + 1.206
               / ((x[nuv2] - 4.62)**2 + 0.263) + Fb)

    AV = RV * ebv
    klamb = RV * (a + b/RV)
    Alamb = AV * (a + b/RV)

    return (Alamb, klamb)


def get_line_map(data_cube, lamb, line_lamb, filter_width=30, cont_width=30,
                 cont_pad=30):
    """
    Calculates a 2D continuum-subtracted emission line map.

    The 2D line map is a collapse of ``data_cube`` along the spectral axis over
    a filter, centred at ``line_mean`` with width ``filter_width``. The
    continuum is estimated by interpolating from two neighbouring wavelength
    windows. The line only map is the filter map minus the interpolated value
    of the continuum.

    Parameters
    ----------
    data_cube : :class:`astropy.io.fits.ImageHDU`
        The data cube extension from which to extract the segmentation map.
    lamb : array-like
        The wavelengths of the datacube spectral axis.
    line_lamb : float
        The central wavelength of the emission line to extract (in units of
        ``data_cube``'s header.
    filter_width : float, optional
        The filter width to extract around the emission line.
    cont_width : float, optional
        The width of the wavelength window to extract either side of the filter
        to define the continuum.
    cont_pad : float, optional
        The padding to apply between the edge of the filter and the start of
        continuum window.

    Returns
    -------
    2d_maps : list
        A length-3 list of 2D maps of the [filter (line+continuum), continuum,
        line only], respectively.
    """
    # Find edges of filter and continuum windows
    low_filt = line_lamb - filter_width/2.
    upp_filt = line_lamb + filter_width/2.
    low_cont1 = low_filt - cont_pad - cont_width
    upp_cont1 = low_filt - cont_pad
    low_cont2 = upp_filt + cont_pad
    upp_cont2 = upp_filt + cont_pad + cont_width
    # Find the nearest sampled wavelengths to our limits
    idx_low_filt = np.abs(lamb - low_filt).argmin()
    idx_upp_filt = np.abs(lamb - upp_filt).argmin() + 1
    idx_low_cont1 = np.abs(lamb - low_cont1).argmin()
    idx_upp_cont1 = np.abs(lamb - upp_cont1).argmin() + 1
    idx_low_cont2 = np.abs(lamb - low_cont2).argmin()
    idx_upp_cont2 = np.abs(lamb - upp_cont2).argmin() + 1
    # Get widths of the windows in order to normalise the
    # continuum to the filter
    filt_width = float(idx_upp_filt - idx_low_filt)
    cont1_width = float(idx_upp_cont1 - idx_low_cont1)
    cont2_width = float(idx_upp_cont2 - idx_low_cont2)
    # Make maps of the three wavelength windows by averaging over the spectral
    # axis
    filter_map = np.nansum(data_cube.data[idx_low_filt:idx_upp_filt, :, :],
                           axis=0)
    cont1_map = np.nansum(data_cube.data[idx_low_cont1:idx_upp_cont1, :, :],
                          axis=0) * (filt_width/cont1_width)
    cont2_map = np.nansum(data_cube.data[idx_low_cont2:idx_upp_cont2, :, :],
                          axis=0) * (filt_width/cont2_width)
    # Determine the mean wavelength of the continuum windows
    cont1_mean = np.average((low_cont1, upp_cont1))
    cont2_mean = np.average((low_cont2, upp_cont2))
    # Interpolate their values to estimate continuum at line mean
    cont_interp = interp1d([cont1_mean, cont2_mean], [cont1_map, cont2_map],
                           axis=0)
    cont_map = cont_interp(line_lamb)
    # Subtract this continuum
    line_map = filter_map - cont_map

    return [filter_map, cont_map, line_map]


def resample_base(basefile, obs_lamb, delta_lamb, buffer_lamb=500):
    """
    Resample starlight base files to match ``delta_lamb`` and restrict coverage
    to limits of ``obs_lamb`` with an additional buffer.
    """
    b = np.genfromtxt(basefile)

    # Extend the ob_lamb array by buffer_lamb on each side, preserving the
    # sampling
    low = obs_lamb[0] - np.arange(1, buffer_lamb/delta_lamb)[::-1] * delta_lamb
    upp = obs_lamb[-1] + np.arange(1, buffer_lamb/delta_lamb) * delta_lamb

    interp_lamb = np.hstack((low, obs_lamb, upp))
    interp_flux = np.interp(interp_lamb, b[:, 0], b[:, 1])
    b_new = np.vstack((interp_lamb, interp_flux)).T
    np.savetxt(basefile, b_new)


def fit_starlight(fargs):
    """
    Perform the fitting of a spectrum with starlight and return the bin number
    and the output staright file.
    """
    bin_num, spec_file, lamb, deltalamb, tmp_dir, bases, low_sn, upp_sn = fargs

    print("fitting bin number {:>5}".format(bin_num), end="\r")
    sys.stdout.flush()
    out_file = os.path.basename(spec_file)+"_out"
    # Write a grid file for starlight (see manual)
    with tempfile.NamedTemporaryFile(prefix="grid_", dir=tmp_dir,
                                     delete=False) as grid_file:
        # we need to specifically set the precision and manually calculate the
        # wavelength window for starlight to model - in some cases the rounding
        # precision losses means our last observed wavelength is not included
        delta_lamb_syn = round(deltalamb, 8)
        lamb_syn_ini = round(lamb[0], 8)
        lamb_syn_fin = round(lamb[-1], 8)+0.01
        # the header
        grid_file.write("\n".join(["1", os.path.join(tmp_dir, bases, ""),
                                   tmp_dir, tmp_dir, tmp_dir, str(RANDSEED),
                                   str(low_sn), str(upp_sn), str(lamb_syn_ini),
                                   str(lamb_syn_fin), str(delta_lamb_syn),
                                   "1.0", "FIT", "1", "1"])+"\n")
        # a line describing this fit containing obs_file config base mask
        # red_law v0 vd out_file
        grid_file.write(" ".join([os.path.basename(spec_file),
                                  "starlight.config", bases+".base",
                                  "starlight.mask", REDLAW, "0.0", "150.0",
                                  out_file])+"\n")

    starlightexe = os.path.join(tmp_dir, "StarlightChains_v04.exe")
    sub = subprocess.Popen("nice {0} < {1} > /dev/null"
                           .format(starlightexe, grid_file.name), shell=True,
                           cwd=tmp_dir)
    sub.wait()
    out_file_path = os.path.join(tmp_dir, out_file)
    try:
        open(out_file_path)
    except IOError:
        out_file = None
    else:
        out_file = out_file_path
    return bin_num, out_file


def fit_emission_lines(fargs):
    bin_num, bin_res, el, vd_init, v0_init, amp_init, stddev_b, off_b, \
        w, filtwidth = fargs
    print("fitting bin number {:>5}".format(bin_num), end="\r")
    sys.stdout.flush()

    spec = bin_res["continuum"]["sl_spec"]
    if w:
        # The starlight output doesn't give uncertainty on regions masked in
        # the fitting (i.e. emission lines!) so we need to copy these over from
        # our original spectrum and account for the normalisation of starlight
        stddev_cube = bin_res["spec"][:, 2]
        stddev_cube[stddev_cube == 0] = np.nan
        spec[:, 3] = bin_res["continuum"]["fobs_norm"]/stddev_cube
    else:
        spec[:, 3] = np.ones(len(spec[:, 0]))

    if filtwidth:
        # filtwidth is the width of a median filter to pass over the residual
        # spectrum to improve the continuum subtraction.
        # This residual function is added to the model spectrum here
        # before fitting emission lines later.
        resid_fn = ndimage.filters.median_filter(spec[:, 1]-spec[:, 2],
                                                 int(round(filtwidth)),
                                                 mode="nearest")
        spec[:, 2] += resid_fn
    else:
        resid_fn = 0.0

    # Contruct a model based on our emission lines
    el_init = _get_emline_model(el).rename(bin_num)
    # Add bounds and ties to the parameters prior to fitting
    for sm in el_init.submodel_names:
        e, wl, n = sm.split("_")
        n = int(n)

        if n == 0:
            # Bound this elem to our offset and stddev limits
            el_init[sm].mean.bounds = [el[e][n] * (1+off_b[0]/ckms),
                                       el[e][n] * (1+off_b[1]/ckms)]
            el_init[sm].stddev.bounds = [el[e][n] * stddev_b[0]/ckms,
                                         el[e][n] * stddev_b[1]/ckms]
        else:
            # If we have a doublet/triplet etc. then tie these to the
            # stddev of the 0th line
            wl0 = el[e][0]
            sm0 = "{}_{:.0f}_0".format(e, wl0)
            el_init[sm].stddev.tied = (lambda x, sm=sm, sm0=sm0:
                                       (x[sm].mean
                                        * (x[sm0].stddev / x[sm0].mean)))
        if sm in ("Halpha_6563_0", "[NII]_6583_0"):
            continue
        # Tie the means of other lines to the anchor forbidden/balmer
        # line as appropriate
        if e.startswith("[") and e.endswith("]"):
            el_init[sm].mean.tied = lambda x, e=e, n=n: (el[e][n]
                                    * x["[NII]_6583_0"].mean/el["[NII]"][0])
        else:
            el_init[sm].mean.tied = lambda x, e=e, n=n: (el[e][n]
                                    * x["Halpha_6563_0"].mean/el["Halpha"][0])

    # Mask regions away from emission lines outside windows of
    # offset_limit + 3*stddev_limit around line rest wavelengths
    rest_lambdas = np.array([_wl for wls in el.values() for _wl in wls])
    lims = 1 + (np.array(off_b) + 3*np.array((-stddev_b[1], stddev_b[1])))/ckms
    fit_limits = rest_lambdas[:, None] * lims
    to_fit = np.zeros(len(spec[:, 0]), bool)
    for lim in fit_limits:
        to_fit[(lim[0] <= spec[:, 0]) & (spec[:, 0] <= lim[1])] = True
    el_spec = spec[to_fit]

    # Use the initial guesses for param combinations later
    if v0_init == "vstar":
        offset_low = 1 + max((bin_res["continuum"]["v0_min"]-200),
                             off_b[0])/ckms
        offset_high = 1 + min((bin_res["continuum"]["v0_min"]+200),
                              off_b[1])/ckms
        offset_init = np.linspace(offset_low, offset_high, 6)
    else:
        offset_init = 1 + np.array(np.array(v0_init)/ckms)
    stddev_init = np.asarray(np.array(vd_init)/ckms)
    # Make combinations of all parameter initial guesses.
    param_comb = list(product(amp_init, offset_init, stddev_init))
    dof = len(fitting._model_to_fit_params(el_init)[0])
    chi2dof = 1e50
    best_fit = None
    # Perform minimisation with LevMar fitter for each combo to find ~global
    # minimum
    for i, comb in enumerate(param_comb, 1):
        # Make an initialisation of the model using this combo
        amp, off, sd = comb
        for sm in el_init.submodel_names:
            elem, wl, n = sm.split("_")
            wavelength = el[elem][int(n)]
            el_init[sm].parameters = [amp, wavelength*off, wavelength*sd]
        # Create the fitter and perform the fit for this combination
        levmar_fitter = fitting.LevMarLSQFitter()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            el_fit = levmar_fitter(el_init, el_spec[:, 0],
                                   el_spec[:, 1]-el_spec[:, 2],
                                   weights=el_spec[:, 3],
                                   maxiter=3000)
        diff = (el_spec[:, 1] - el_spec[:, 2]) - el_fit(el_spec[:, 0])
        _chi2 = np.nansum((diff * el_spec[:, 3])**2)
        _chi2dof = _chi2/dof
        if _chi2dof < chi2dof:
            chi2dof = _chi2dof
            best_fit = el_fit
            # Here we explicitly pass the fit info dict to our fitted model
            best_fit.fit_info = levmar_fitter.fit_info
            best_fit.chi2dof = chi2dof

    # We cannot put constraints on the parameters and still get a covariance
    # matrix, so here we correct any negative amplitudes since these are
    # suppose to be emission features
    for sm in best_fit.submodel_names:
        if best_fit[sm].amplitude.value <= 0:
            best_fit[sm].amplitude.value = 0

    # Save all parameters for the model and uncertainties.
    res = _model_to_res_dict(best_fit, el, bin_res["continuum"]["fobs_norm"],
                             filtwidth, resid_fn)

    return bin_num, res


def parse_starlight(sl_output):
    """
    Returns a dictionary of useful quantities from a starlight output file.
    """
    results_dict = {}
    try:
        lines = open(sl_output).readlines()
    except IOError:
        print("cannot open {}!".format(sl_output))
        return None
    n_lines = len(lines)
    if n_lines < 20:
        return None
    spec_line = None
    for i, line in enumerate(lines):
        if not line.strip():
            continue
        m = re.search(r"\[(.*?)\]", line)
        if m:
            val = m.group(1).split()[0]
            if val in ["fobs_norm", "Mini_tot", "Mcor_tot", "v0_min",
                       "vd_min", "AV_min", "chi2/Nl_eff", "adev",
                       "NOl_eff", "Nl_eff", "N_base", "Nl_obs"]:
                results_dict[val] = float(line.split()[0])
                if val == "Nl_obs":
                    spec_line = i + 1
        if line.startswith("# j     x_j(%)"):
            base_line = i + 1
    if spec_line is None:
        return None
    results_dict["sl_spec"] = np.genfromtxt(sl_output, skip_header=spec_line)
    results_dict["bases"] = np.genfromtxt(sl_output, skip_header=base_line,
                                          max_rows=results_dict["N_base"],
                                          usecols=(0,1,2,3,4,5,6,8))
    return results_dict


def _model_to_res_dict(model, el, fobs_norm, filtwidth, resid_fn):
    """
    Given an emission line fitted model, return a dictionary of fitted
    parameters and uncertainties.
    """
    res = {}
    res["chi2dof"] = model.chi2dof
    res["filtwidth"] = filtwidth
    res["resid_fn"] = resid_fn * fobs_norm
    try:
        fitted_uncerts = np.diag(model.fit_info["param_cov"])**0.5
    except ValueError:
        print("no covariance matrix computed for bin {}, cannot"
              " compute fit uncertainties".format(model.name))
        dof = len(fitting._model_to_fit_params(model)[0])
        fitted_uncerts = np.full(dof, np.nan)
        res["bad"] = 1
    else:
        res["bad"] = 0

    # Find the location of our balmer and forbidden offset anchor
    # lines in the uncert array
    param_idx = np.array(fitting._model_to_fit_params(model)[1])
    sm_balmer_idx = model.submodel_names.index("Halpha_6563_0")
    sm_forbidden_idx = model.submodel_names.index("[NII]_6583_0")
    uncertb_idx = np.argwhere(param_idx == sm_balmer_idx*3)
    uncertf_idx = np.argwhere(param_idx == sm_forbidden_idx*3)

    res["lines"] = {}
    j = 0
    for i, sm in enumerate(model.submodel_names):
        e, wl, n = sm.split("_")
        n = int(n)
        name = "{}_{}".format(e, wl)  # `elem_wl`
        res["lines"][name] = {}
        sm_res = res["lines"][name]
        # Store the amplitude, mean and stddev of the gaussians
        sm_res["fit_params"] = model.parameters[i*3:i*3+3]
        sm_res["fit_params"][0] *= fobs_norm  # remove amp normalisation
        # Retrieve uncertainties using indexes of our
        # primary balmer/forbidden lines since the fitting
        # does not return uncertainties for tied parameters
        sm_uncerts = np.empty(3)
        sm_uncerts[0] = fitted_uncerts[j] * fobs_norm  # amplitude
        if sm in ("Halpha_6563_0", "[NII]_6583_0"):
            j += 1
            sm_uncerts[1] = fitted_uncerts[j]  # mean
            j += 1
            sm_uncerts[2] = fitted_uncerts[j]  # stddev
        else:
            if e.startswith("[") and e.endswith("]"):
                sm_uncerts[1] = fitted_uncerts[uncertf_idx+1]  # mean
            else:
                sm_uncerts[1] = fitted_uncerts[uncertb_idx+1]  # mean
            if n == 0:
                j += 1
                sm_uncerts[2] = fitted_uncerts[j]  # stddev
            else:
                if sm == "[NII]_6583_1":
                    sm_uncerts[2] = fitted_uncerts[uncertf_idx+2]  # stddev
                else:
                    wl0 = el[e][0]
                    sm0 = "{}_{:.0f}_0".format(e, wl0)
                    sm0_idx = model.submodel_names.index(sm0)
                    uncert0_idx = np.argwhere(param_idx == sm0_idx*3)
                    sm_uncerts[2] = fitted_uncerts[uncert0_idx+1]  # stddev

        # Store the rest wavelength to calculate offsets later
        sm_res["rest_lambda"] = el[e][n]
        sm_res["fit_uncerts"] = sm_uncerts
        j += 1

    return res


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name="shiftedcmap"):
    """
    Function to offset the "center" of a colormap. Useful for data
    with a negative min and positive max and you want the middle of
    the colormap's dynamic range to be at zero

    Parameters
    ----------
    cmap : matplotlib.cm
        The matplotlib colormap to be altered
    start : float
        Offset from lowest point in the colormap's range. Should be between 0.0
        and ``midpoint``.
    midpoint : float
        The new center of the colormap. Should be between 0.0 and 1.0. In
        general, this should be 1 - vmax/(vmax + abs(vmin)) For example if your
        data range from -15.0 to +5.0 and you want the center of the colormap
        at 0.0, ``midpoint`` should be set to 1 - 5/(5 + 15)) or 0.75
    stop : float
        Offset from highest point in the colormap's range. Should be between
        ``midpoint`` and 1.0.

    """
    cdict = {
        "red": [],
        "green": [],
        "blue": [],
        "alpha": []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict["red"].append((si, r, r))
        cdict["green"].append((si, g, g))
        cdict["blue"].append((si, b, b))
        cdict["alpha"].append((si, a, a))

    newcmap = colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


if __name__ == "__main__":
    pass
