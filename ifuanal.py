"""
IFUANAL

For the analysis of IFU data cubes.
"""

__version__ = "0.1.0"
__author__ = "J. Lyman"

from itertools import repeat, cycle, product
import json
from math import log10, pi
import os
import re
import tempfile
import shutil
import subprocess
import sys
import warnings

from astropy.constants import c
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.modeling import models, fitting
import astropy.units as u
import astropy.wcs as wcs
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import gridspec, ticker, cm, colors, rc
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
from scipy.interpolate import interp1d

from voronoi import voronoi
import dill as pickle
import pathos.multiprocessing as mp

# The file and directory path to this script
FILEPATH = os.path.realpath(__file__)
FILEDIR = os.path.dirname(FILEPATH)

# Colours for plots
OBSCOL = "k" # observed spectrum
SYNTHCOL = "#D95F02" # starlight fit
MASKCOL = "#7570B3" # masked regions in fit
CLIPPEDCOL = "red" # clipped by starlight
ERRCOL = "#666666" # stddev
ZCOL = np.array(["#1B9E77", "#E7298A", "#a6761d","#66A61E", "#E6AB02", 
                 "#7570B3"]) # metallicities

# Random seed for starlight grid file
RANDSEED = 999
# Reddening law to use in starlight
REDLAW = "CCM"

# Speed of light in km/s
ckms = c.to("km/s").value 

# Make plots prettier
rc('font',**{'family':'serif','serif':['Times New Roman'],'size':14})
rc('text', usetex=True)
rc('image', cmap='viridis')

#TODOs

#TODO    : generalise the emission line model to any emission_line.json input.
#          add numbering to components in emission_lines.json with Halpha 0 and
#          NII 1, then split by balmer and forbidden?  When constructing
#          emissionline model need to rename them to the name of the lines?

#TODO    : bin_inspection - i.e.: make a plot of all indiv specs in bin and the
#          weighted spec, show location on voronoi plot, print results of
#          fitting, overplot fit on weighted spec etc.

#TODO    : make all bin numbers, the x/y coordinates, the weighted spec, the
#          starlight outputs the error spec, the s/n of the bin etc. all in a
#          big astropy table. Then can parse this easily for output analysis

#TODO    : put all output into a subdir of the basedir instead

#TODO    : output plots with RA/DEC instead of just pixels?

#TODO    : add an option to the custom bin method that will simulate an SDSS
#          fibre on the host.

#TODO    : only calculate metallicities when its within bound of validity for
#          the indicator (otherwise nan), and also propogate errors onto the
#          value (excluding intrinsic scatter about the relations themselves).

#TODO:     how to handle emission models that are fit but on limits of bounds
#          correctly?  this produces a model nan uncerts. Should these be
#          thrown out? Some testing needed to see if bounds are cutting out
#          crap while not allowing good to make stupid fits.

#TODO:     expose emission_lines as an argument to IFUCube so user can give own
#          emission_lines.json file

#TODO:     create data/metallicity_indicators.json for custom indicators to be
#          used/added

class IFUCube(object):
    """
    Parent to be called by initialisation of child.

    Parameters
    ----------
    data_cube : :class:`astropy.io.fits.ImageHDU`
        The data extension of the IFU FITS file.
    stddev_cube : :class:`astropy.io.fits.ImageHDU`
        The uncertainty (standard deviation) extension of the IFU FITS file.
    base_name : str
        The base name to be used for all output.
    RV : float, optional
        The RV value to use for the CCM extinction law.
    vor_sn : float, optional
        The target signal-to-noise of the voronoi binning algorithm.
    sl_dir : None or str, optional
        The directory containing starlight files and bases. The default
        ``None`` will use the `starlight/` subdirectory.
    """
    def __init__(self, data_cube, stddev_cube, base_name, RV=3.1, vor_sn=20,
                 sl_dir=None):

        self.data_cube = data_cube
        self.stddev_cube = stddev_cube
        self.base_name = base_name

        self.RV = RV
        self.data_shape = data_cube.data.shape
        self.stddev_shape = stddev_cube.shape
        if self.data_shape != self.stddev_shape:
            raise ValueError("data_cube and stddev_cube must have same shape!")
        
        # Make a wavelength array and store its units
        sl = data_cube.header["CRVAL3"] # start lambda
        rl = data_cube.header["CRPIX3"] # ref pix lambda
        dl = data_cube.header["CD3_3"] # delta lambda
        self.lamb = sl + (np.arange(self.data_shape[0]) - (rl - 1)) * dl
        self.delta_lamb = dl
        self.lamb_units = u.Unit(data_cube.header["CUNIT3"])

        self.vor_sn = vor_sn # desired S/N of voronoi bins

        if sl_dir is None:
            self.sl_dir = os.path.join(FILEDIR, "starlight")
        elif not os.path.isdir(sl_dir):
            raise AttributeError("{} is not an accesible "
                                 "directory".format(sl_dir))
        else:
            self.sl_dir = sl_dir
        with open(os.path.join(FILEDIR, "data", "emission_lines.json")) as f:
            self.emission_lines = json.load(f)

        self.nucleus = None
        self.n_cpu = int(min(mp.cpu_count()-1,mp.cpu_count()*0.9))

        self.sl_output = {} # dictionary of {[bin number]: ["starlight
                            # outfile", "spec infile",}
        self.emline_output = {} # dictionary of {[bin number]: [emline_model,
                                # emline_model chi2dof, good/bad fit]}
        self.results = None # This will hold everything you could possibly want

        self.emline_model = None # The compound model to use in fitting
                                 # emission lines

    def deredden(self):
        """
        Correct data and stddev cubes for Galactic dust.

        Based on CCM law and the value of the header card `IFU_EBV`.
        After correction, the header card `IFU_EBV` is set to zero.

        """
        ebv = self.data_cube.header["IFU_EBV"]
        if ebv == 0:
            print "ebv = 0, skipping deredden()"
            return
        print("dereddening with E(B-V) = {:.3f}mag and RV = {}"
              .format(ebv, self.RV))
        Alamb = get_Alamb(self.lamb, ebv, self.RV, self.lamb_units)
        corr = 10**(0.4 * Alamb)

        # Multiply our flux and stddevs by the correction
        self.data_cube.data *= corr[:, None, None]
        self.stddev_cube.data *= corr[:, None, None]
        # The cube is now dereddened
        self.data_cube.header["IFU_EBV"] = 0.

    def deredshift(self):
        """
        Correct wavelengths of the data and stddev cubes to restframe.

        After correction, the header card `IFU_Z` is set to zero.
        """
        z = self.data_cube.header["IFU_Z"]
        if z == 0:
            print("redshift = 0, skipping deredshift()")
        print("deredshifting from z = {}".format(z))
        # Update the lambda array
        self.lamb /= 1. + z
        self.delta_lamb /= 1 + z
        # Update data header
        self.data_cube.header["CRVAL3"] /= 1. + z # start lambda
        self.data_cube.header["CRPIX3"] /= 1. + z # ref pix lambda
        self.data_cube.header["CD3_3"]  /= 1. + z # delta lambda
        # Update stddev header
        self.stddev_cube.header["CRVAL3"] /= 1. + z 
        self.stddev_cube.header["CRPIX3"] /= 1. + z 
        self.stddev_cube.header["CD3_3"]  /= 1. + z
        # The cube is now restframe
        self.data_cube.header["IFU_Z"] = 0.

    def mask_regions(self, centres, r):
        """
        Set spaxels in circular regions defined by ``centres`` with radius
        ``r`` to NaNs.

        Parameters
        ----------
        centres : list of lists
            Coordinates of mask centres of the form ((x1,y1),...(xn,yn)).
            Should be given as zero-indexed values.
        r : int or array_like
            The radius (radii) of the masks. If given as ``int``, same radius
            used for all masks given in ``centres``. If ``array_like`` then this
            will be cycled over for all ``centres``.
        """
        print("masking regions")
        if isinstance(r, int):
            r = [r]
        elif len(r) != len(centres) and len(r) != 1:
            warnings.warn("Number of radii ({}} not same as centres ({}). "\
                          "Radii will be cycled over.", RuntimeWarning)
        for centre, _r in zip(centres, cycle(r)):
            x_cen, y_cen = centre
            y, x = np.ogrid[-y_cen:self.data_shape[1]-y_cen,
                            -x_cen:self.data_shape[2]-x_cen]
            mask = x*x + y*y <= _r*_r
            self.data_cube.data[:,mask] = np.nan

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
            a gaussian in pixels. If 0 then will for nucleus as ``xc``, ``yc``.
        lamb_low, lamb_upp : float, optional
            The wavelength limits over which to sum.
        plot : bool, optional
            Make a plot of the data, model and residual).
        """

        if usewcs:
            w = wcs.WCS(self.data_cube.header)
            s = SkyCoord(xc, yc, unit=(u.hourangle, u.deg))
            xc, yc = s.to_pixel(w)
            if np.any(np.isnan([xc,yc])):
                raise AttributeError("Couldn't find pixel location for {}"\
                                     .format(s.to_string(hmsdms)))

        if box_size == 0:
            self.nucleus = (xc, yc)
            print("set nucleus as {}".format(self.nucleus))
            return
        elif not box_size <= xc <= self.data_shape[2] - box_size \
             or not box_size <= yc <= self.data_shape[1] - box_size:
            raise AttributeError("box must be fully within the image, use "
                                 "box_size=0 to force a location outside "
                                 "the FOV.")
            return 

        # Find the nearest sampled wavelengths to our limits
        idx_low = np.abs(self.lamb - lamb_low).argmin()
        idx_upp = np.abs(self.lamb - lamb_upp).argmin() + 1

        # Cutout the subarray then sum and normalise it
        cutout = self.data_cube.data[:,yc-box_size:yc+box_size+1,
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

        self.nucleus = (round(xc-box_size+g.x_mean,3), 
                        round(yc-box_size+g.y_mean,3))
        print("found nucleus as {}".format(self.nucleus))

        # Plot the data with the best-fit model
        if plot:
            plt.close("all")
            plt.figure(figsize=(9.5, 3))
            for i, vt in enumerate([(z, "datacube"), 
                                    (g(x,y), "model"), 
                                    (z-g(x,y), "residual")], 1):
                vals, title = vt
                plt.subplot(1, 3, i)
                plt.imshow(vals, origin='lower', interpolation='nearest', 
                           vmin=np.min(z), vmax=1.)
                plt.autoscale(False)
                plt.plot(g.x_mean, g.y_mean, "kx", markersize=10)
                plt.gca().get_xaxis().set_visible(False)
                plt.gca().get_yaxis().set_visible(False)
                plt.title(title)
            plt.gcf().text(0.5,0.05,"nucleus = "+str(self.nucleus), 
                           ha="center", va="center")
            plt.savefig(self.base_name+"_nucleus.pdf", bbox_inches="tight")
       

    def voronoi_bin(self, lamb_low=5590., lamb_upp=5680., cont_lamb_low=None,
                    cont_lamb_upp=None, clobber=False, min_sn=3):
        """
        Apply the voronoi binning algorithm to the data cube.

        The target signal-to-noise (S/N) is given by the class attribute
        ``vor_sn``. S/N calculation is performed in the wavelength window given
        by ``lamb_low`` and ``lamb_upp``. Only spaxels that have S/N >=
        ``min_sn`` will be considered for binning. This will write results in a
        file with the suffix `_voronoibins.txt`.

        Emission line signal to noise can be calculated by specifying
        ``cont_lamb_low`` and ``cont_lamb_upp``. In this case the S/N in this
        window is removed from the ``lamb_low`` and ``lamb_upp`` window to give
        the S/N in the emission lines.

        Parameters
        ----------
        lamb_low, lamb_upp : float, optional
            The wavelength limits over which to calculate the S/N.
        cont_lamb_low, cont_lamb_upp : float or None optional
            The wavelength limits of the continuum to subtract in calculating
            S/N (defaults to ``None``, ``None`` - i.e. no subtraction)
        clobber: bool, optional
            Whether to overwrite an existing output. Searches for the file
            ``base_name``+`_voronoibins.txt`, if this exists and ``clobber`` = 
            ``False`` then just use that.
        min_sn : float, optional
            The minimum S/N of a spaxel to be considered for the algorithm.

        Example
        -------
        To produce bins determined by their signal to noise in H\
        :math:`\\alpha` + [NII]:

        ``>>> IFUCube.voronoi_bin(lamb_low=6540., lamb_upp=6580.,``\
        ``cont_lamb_low=6600, cont_lamb_upp=6640)``


        References
        ----------
        Uses the algorithm of [CC03]_.

        .. [CC03] Cappellari, M & Copin, Y, "Adaptive spatial binning of 
           integral-field spectroscopic data using Voronoi tessellations",
           MNRAS, 2003.

        """
        print("voronoi_binning with S/N target of {}".format(self.vor_sn))
        vor_output_file = self.base_name+"_voronoibins.txt"

        if clobber is False:
            # Search for an existing output
            try:
                vor_output = np.genfromtxt(vor_output_file)
            except IOError, ValueError:
                print("no valid voronoi binning file found")
            else:
                # Set attributes based on the existing file
                self.bin_nums = np.sort(np.unique(vor_output[:,3]))\
                                  .astype("int")
                self.x_bar, self.y_bar = vor_output[:,5], vor_output[:,6]
                self.bin_sn = vor_output[:,4]
                self.vor_output = vor_output
                return

        if not all((self.lamb[0] < lamb_low < self.lamb[-1], 
                    self.lamb[0] < lamb_upp < self.lamb[-1])):
            raise ValueError("S/N window not within wavelength range of cube")
        # Find the nearest sampled wavelengths to our limits
        idx_low = np.abs(self.lamb - lamb_low).argmin()
        idx_upp = np.abs(self.lamb - lamb_upp).argmin() + 1
        
        # Sum the data and stddev cubes to get signal and noise while
        # catching the warning when spaxels only have nans
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
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
        vor_input = vor_input[~np.isnan(vor_input[:,3])]
        #  remove any with negative or zero noise
        vor_input = vor_input[vor_input[:,3] > 0]
        #  also reduce to only spaxels with S/N >= min_sn
        vor_input = vor_input[(vor_input[:,2]/vor_input[:,3]) >= min_sn]

        # Split columns as voronoi wants separate inputs
        x, y, sig, noi = vor_input[:, [0,1,2,3]].T

        # Call the voronoi binning script
        vor_plot = self.base_name + "_voronoibins.pdf"
        res  = voronoi.voronoi_2d_binning(x, y, sig, noi, targetSN=self.vor_sn,
                                          cvt=True, pixelsize=1, plot=vor_plot,
                                          quiet=False, n_cpu=self.n_cpu)
        bin_num, x_node, y_node, x_bar, y_bar, bin_sn, n_pix, scale = res

        # Save the output of voronoi binning
        vor_output = np.column_stack([x, y, sig/noi, bin_num, bin_sn[bin_num],
                                      x_bar[bin_num], y_bar[bin_num]])
        with open(vor_output_file, "wb") as f:
            f.write(b"#  x    y  sn_spec bin_num   sn_bin bin_xbar bin_ybar\n")
            np.savetxt(f, vor_output, 
                       fmt=b"%4i %4i %8.3f %7i %8.3f %8.3f %8.3f")

        # Some class attributes to assign
        self.bin_nums = np.sort(np.unique(bin_num)).astype("int")
        self.x_bar, self.y_bar = x_bar, y_bar
        self.bin_sn = bin_sn
        self.vor_output = vor_output


    def emission_line_bin(self, min_peak_flux, min_frac_flux, max_radius,
                          min_flux, **kwargs):
        """
        Apply the HII explorer [SFS]_ binning algorithm to the datacube.

        This method will bin spaxels by attempting to determine distinct
        emission line regions. An emission line map (usually H:math:`\\alpha`)
        is created by :func:`~ifuanal.get_line_map`` through subtraction of a
        continuum from a narrow band filter centred on the emission line. Peaks
        above ``min_peak_flux`` in the emission line map are seeds for bins
        which are expanded so long as neighbouring spaxels are above
        ``min_frac_flux`` of the peak and within ``max_radius``
        distance. Pixels below ``min_flux`` are excluded from binning.

        See :func:`~ifuanal.get_line_map` for more information on the kwargs
        used to define the wavelength windows of the line and continuum.

        Parameters
        ----------
        min_peak_flux : float
            The minimum flux (in what units?!) for a spaxel to be considered as
            a new bin seed.
        min_frac_flux : float
            The minimum flux of a spaxel, as a fraction of the bin's peak flux,
            to be considered a member of the bin.
        max_radius : float
            The maximum radius (in what units?!) allowed for a bin.
        min_flux : float
            The minimum flux of a spaxel for consideration in the binning
            routine.
        line_mean : float
            The central wavelength of the emission line to extract (in units of
            ``data_cube``'s header.
        filter_width : float
            The filter width to extract around the emission line.
        cont_width : float
            The width of the wavelength window to extract either side of the
            filter to define the continuum.
        cont_pad :
            The padding to apply between the edge of the filter and the start
            of continuum window.

        References
        ----------
        .. [SFS] S.F. Sanchez, HII_explorer,
           http://www.caha.es/sanchez/HII_explorer/
        """
        pass #TODO #TODO #TODO #TODO #TODO 

    def add_custom_bin(self, centre, r):

        """
        Create a custom bin to analyse in addition to the vonoroi bins.

        Custom bins have negative values (beginning at -1) in all output etc.

        Parameters
        ----------
        centre : array_like
            Length 2 array giving the x,y centre of the bin
        r : int
            The radius of the bin in pixels.


        Example
        -------
        Creating an SDSS-like  3 arcsec fibre centred on the host:

        1. set the nucleus with :meth:`set_nucleus()` (it must be in the FOV)
        2. determine the ``plate_scale`` of the cube (arcsec/pixel) 
        3. ``>>> cube.add_custom_bin(cube.nucleus, 3/plate_scale)``
        """

        x_cen, y_cen = centre 
        y, x = np.ogrid[-y_cen:self.data_shape[1]-y_cen, 
                        -x_cen:self.data_shape[2]-x_cen]
        dataidx = np.where(x*x + y*y <= r*r)

        voridx = np.zeros(len(dataidx[0]), dtype="int")
        for i,di in enumerate(zip(dataidx[0], dataidx[1])):
            voridx[i] = np.where((self.vor_output[:,:2] == di)\
                                 .all(axis=1))[0][0]

        # Make a new row of data that we will append to the existing
        # output from the voronoi binning
        vor_append = self.vor_output[voridx]
        vor_append[:,4] = -999
        bin_num = -1
        while True:
            if not np.any(self.vor_output[:,3] == bin_num):
                break
            bin_num -= 1
            continue
        vor_append[:,3] = bin_num
        self.bin_nums = np.append(self.bin_nums, bin_num)
        density = vor_append[:,2]**2
        mass = np.sum(density)
        vor_append[:,5] = np.sum(vor_append[:,0]*density)/mass
        vor_append[:,6] = np.sum(vor_append[:,1]*density)/mass     

        self.vor_output = np.vstack((self.vor_output, vor_append))
        print("added custom bin {} to the list".format(bin_num))

    def get_bin_num(self, loc):
        """
        Return the bin number for a given location.

        Parameters
        ----------
        loc : array_like
            Length 2 array giving the (``x``, ``y``) to get bin number for.

        Returns
        -------
        bin_num : int
            The bin number that contains ``loc``. Returns None if not in any bin.
        """
        x, y = map(int, loc)
        idx = np.where((self.vor_output[:,:2] == [x,y]).all(axis=1))[0]
        if idx.size == 1:
            return np.asscalar(self.vor_output[idx,3])
        else:
            print("didn't find a bin at location {}, {}".format(x,y))
            return None


    def get_single_spectrum(self, x, y):
        """
        Return the spectrum for spaxel at location ``x``, ``y``.

        Returned array is shape (N,4), where N is number of wavelength
        samples of cube and columns are wavelength, flux, flux_err, flag.
        Flag = 2 for negative or nan fluxes (dummy values for flux and flux_err
        are inserted for these wavelengths).
        """
        spec = np.empty((self.lamb.size, 4))
        spec[:,0] = self.lamb

        spec[:,1] = self.data_cube.data[:,y,x] # numpy axes switch
        spec[:,2] = self.stddev_cube.data[:,y,x]
        # STARLIGHT ignores flags >=2
        bad_idx = np.isnan(spec[:,1]) | np.isnan(spec[:,2])
        spec[:,3] = (bad_idx).astype("int") * 2

        return spec

    def get_weighted_spectrum(self, x, y):
        """
        Return the weighted median spectrum for spaxels at locations ``x``, ``y``.

        Similar to ``get_single_spectrum`` except ``x`` and ``y`` are arrays.
        The single spectra given by these locations are combined using the
        weighted mean of the fluxes. Returns array of same form as 
        ``get_single_spectrum``.
        """
        # Use weighted arthimetic mean and variance of weighted mean
        # arrays to hold all the flux and flux_stddev values in spaxels
        spaxels_flux = self.data_cube.data[:,y,x] # numpy axes switch
        spaxels_stddev = self.stddev_cube.data[:,y,x]

        # Find any wavelengths where we have >75% spaxels as nans 
        # and flag that wavelength as bad for starlight
        bad_idx = np.isnan(spaxels_flux) | np.isnan(spaxels_stddev)
        num_bad = np.sum(bad_idx, axis=1)
        bad_lamb = num_bad > 0.75 * x.size

        # Array to hold final weighted-mean spectrum - same format 
        # as get_single_spectrum() 
        spec = np.empty((self.lamb.size, 4))
        spec[:,0] = self.lamb
        # Use masked arrays to cover the nans while preserving shape
        spaxels_flux_ma = np.ma.masked_array(spaxels_flux, bad_idx)
        spaxels_stddev_ma = np.ma.masked_array(spaxels_stddev, bad_idx)
        # Calculate the weighted mean
        spec[:,1] = np.ma.average(spaxels_flux_ma,
                                  weights=1/spaxels_stddev_ma**2, axis=1)
        # Calculate the weighted stddev
        diff = spaxels_flux_ma - spec[:,1][:,None]
        var = np.ma.average(diff**2, weights=1/spaxels_stddev_ma**2, axis=1)
        spec[:,2] = np.ma.sqrt(var)
        # STARLIGHT ignores flags >=2
        spec[:,3] = bad_lamb.astype("int") * 2

        return spec

    def yield_spectra(self, bin_num):
        """
        Yield the spectrum for a given ``bin_num``.

        If ``bin_num`` is an array, produce a generator to obtain spectra of
        each. Calls ``get_single_spectrum`` or ``get_weighted_spectrum`` as
        appropriate (i.e. if number of spaxels in bin is 1 or >1)
        """
        try:
            self.vor_output
        except AttributeError:
            raise("No voronoi output found. Run voronoi_bin() first")

        n_bins = len(bin_num)
        for i,bn in enumerate(bin_num,1):
            # Get the indices and number of spaxels that are in our chosen bin
            spaxels_idx = self.vor_output[:,3] == bn
            spaxels_num = np.sum(spaxels_idx)
            print("Calculating weighted mean of {:>3} spaxels in bin {:>5} "
                  "({:>5}/{:>5})\r".format(spaxels_num, bn, i, n_bins)),
            sys.stdout.flush()

            if spaxels_num == 0:
                raise ValueError("No spaxels found in bin number: "
                                 "{}".format(bin_num))
            elif spaxels_num == 1:
                # Just return the single spaxel
                x, y = self.vor_output[spaxels_idx,:2][0]
                spec = self.get_single_spectrum(x, y)
            else:
                # Get the x and y positions of our spaxels in the data cube
                x = self.vor_output[spaxels_idx,0].astype("int")
                y = self.vor_output[spaxels_idx,1].astype("int")
                # and get the weighted spectrum of these positions
                spec = self.get_weighted_spectrum(x,y)        
            yield spec



    def run_starlight(self, bin_num=None, base_name="bc03", lamb_low=5590., 
                      lamb_upp=5680., use_tmp_dir=None, clobber=False,
                      append=False):
        """
        Run the starlight fitting of spectra in the cube.
        
        This method will create a file with the suffix `_sloutputs.txt` which
        gives the bin numbers and their corresponding starlight output files.
        If a run has been previously performed, the contents of this file can
        be read - however the files containing the starlight output must
        also exist.
        
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
            must be within the wavelengh range of the cube anyway (defaults 
            to 5590.->5680.).
        use_tmp_dir: None or str, optional
            The temporary directory where the resampled bases and output will
            be stored. If None a safe temporary directory will be made. A
            previous temporary directory used can be given and the resampled
            bases in there will be used (defaults to None).
        clobber : bool, optional
            Whether to overwrite the `_sloutput.txt` output if it exists. If
            this file exists and ``clobber = False`` then its contents are read
            and ``parse_starlight`` will use the pre-existing output files
            (defaults to False).
        append : bool, optional
            NOT IMPLEMENTED YET. (Allow method to overwrite sl_output only for
            ``bin_num`` and append these to the output rather than overwrite) -
            NOTE sl_output is set as empty dict on __init__)
        """
        print("running starlight fitting")

        if not clobber and not append:
            try:
                f = open(self.base_name+"_sloutputs.txt", "r")
            except IOError:
                print("didn't find existing starlight output file.")
            else:
                print("found {}".format(self.base_name+"_sloutputs.txt"))
                self.sl_output = {}
                for line in f:
                    if line.startswith("#"):
                        continue
                    bn, of, sf = line.split()
                    if of == "None":
                        of = None
                    self.sl_output[int(bn)] = [of, sf]
                return
        elif clobber and append:
            print("can't `clobber` *and* `append`.")
            return
        elif not clobber and append:
            pass
            #TODO

        if bin_num is None:
            bin_num = self.bin_nums
        elif bin_num not in self.bin_nums:
            ValueError("Bin number {} does not exist".format(bin_num))
        elif isinstance(bin_num, (int,float)):
            bin_num = [bin_num]
        print("fitting {} bins...".format(len(bin_num)))
        if use_tmp_dir:
            sl_tmp_dir = os.path.join(use_tmp_dir,"")
            print("using STARLIGHT tmp directory {}".format(sl_tmp_dir))
        else:
            # We need to copy everything to a temporary folder since
            # starlight has filepath character length issues, plus 
            # then it doesn't clog up cwd.
            # Make a safe temporary directory to hold everything
            sl_tmp_dir = os.path.join(tempfile.mkdtemp(prefix="starlight_"),"")
            print("STARLIGHT tmp directory for this run is {}"
                  .format(sl_tmp_dir))
            # Copy over the bases and resample them to that of our data
            bases_tmp_dir = os.path.join(sl_tmp_dir, base_name)
            # Because shutil.copytree requires `dst` doesn't exist:
            shutil.rmtree(sl_tmp_dir) 
            shutil.copytree(self.sl_dir, sl_tmp_dir)
            d = os.path.join(sl_tmp_dir, base_name)
            basefiles = [os.path.join(d, f) for f in os.listdir(d)]
            Nbasefiles = len(basefiles)
            for i, basefile in enumerate(basefiles,1):
                print("resampling base files {}/{}\r".format(i, Nbasefiles)),
                sys.stdout.flush()
                resample_base(basefile, self.lamb, self.delta_lamb)
        print
        # Compute the spectra and write the spectrum to a temp file
        spectra = self.yield_spectra(bin_num)
        spec_files = []
        for spec in spectra:
            with tempfile.NamedTemporaryFile(prefix="spec_",dir=sl_tmp_dir,
                                             delete=False) as spec_file:
                np.savetxt(spec_file, spec, fmt="%14.8f")
                spec_files.append(spec_file.name)
        print
        # multiprocessing params for starlight pool of workers
        p = mp.Pool(self.n_cpu)
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
        print("Starlight fits of {} bins complete".format(len(bin_num)))

        #TODO
        #if append:
        #    file_mode = "a"
        #else:
        #    file_mode = "w"
        #if clobber:
        #    self.sl_output = {}
        with open(self.base_name+"_sloutputs.txt", "w") as f: #"w"->file_mode?
            f.write("# bin starlight_file\n")
            for bn, of, sf in bin_out_files:
                self.sl_output[int(bn)] = [of, sf]
                f.write("{:>5} {} {}\n".format(bn, of, sf))


    def parse_results(self, clobber=False):
        """
        Creates and pickles a dictionary containing all results.

        The ``results`` dictionary is created by parsing the starlight output
        files (see ``run_starlight``) and copying other properties (e.g.
        ``nucleus``). The dictionary is pickled to file with suffix 
        `_results.pkl`.

        Parameters
        ----------
        clobber : bool, optional
            If ``clobber`` is not ``True`` and ``results`` exists, then it will not
            do the parsing.
        """

        print("parsing results")
        if self.results is not None and not clobber:
            print("`results` dict exists, use clobber=True to overwrite")
            return

        try:
            self.vor_output
        except AttributeError:
            raise("no voronoi output found. run voronoi_bin() first")


        self.results = {}
        self.results["nucelus"] = self.nucleus

        # Populate all bin entries in advance
        self.results["bin"] = dict.fromkeys(self.bin_nums)
        n_bins = len(self.bin_nums)
        for i, bn in enumerate(self.bin_nums, 1):
            print("parsing starlight output {}/{}\r".format(i, n_bins)),
            sys.stdout.flush()
            spaxels_idx = self.vor_output[:,3] == bn
            if self.sl_output[bn][0] is not None:
                self.results["bin"][bn] = parse_starlight(self.sl_output[bn][0])
                if self.results["bin"][bn] is None: #i.e. failed to get Nl_Obs
                    print("failed to find [Nl_Obs] line, or suspiciously "
                          "few lines in {}".format(self.sl_output[bn][0]))
                    self.results["bin"].pop(bn)
                    self.bin_nums = np.delete(self.bin_nums,
                                              np.argwhere(self.bin_nums == bn))
                    continue
                self.results["bin"][bn]["spec"] = np.genfromtxt(
                                                      self.sl_output[bn][1])
                x_spax = self.vor_output[spaxels_idx,0].astype("int")
                y_spax = self.vor_output[spaxels_idx,1].astype("int")
                self.results["bin"][bn]["x_spax"] = x_spax
                self.results["bin"][bn]["y_spax"] = y_spax
                x_bar = self.vor_output[spaxels_idx,5][0]
                y_bar = self.vor_output[spaxels_idx,6][0]
                self.results["bin"][bn]["x_bar"] = x_bar
                self.results["bin"][bn]["y_bar"] = y_bar
                if self.nucleus is not None:
                    nx, ny = self.nucleus[0], self.nucleus[1]
                    distances = ((x_spax - nx)**2 
                                 + (y_spax - ny)**2)**0.5
                    self.results["bin"][bn]["dist_min"] = np.min(distances)
                    self.results["bin"][bn]["dist_max"] = np.max(distances)
                    self.results["bin"][bn]["dist_bar"] = ((x_bar - nx)**2 +
                                                           (y_bar - ny)**2)**0.5

                else:
                    self.results["bin"][bn]["distance_min"] = None
                    self.results["bin"][bn]["distance_max"] = None
                    self.results["bin"][bn]["distance_bar"] = None
            else:
                print("no valid sl_output found for {}, removing".format(bn))
                self.results["bin"].pop(bn)
                self.bin_nums = np.delete(self.bin_nums,
                                          np.argwhere(self.bin_nums == bn))     

    def plot_starlight_results(self, bin_num):
        """
        Plot the spectrum and results of starlight continuum fitting for
        ``bin_num``

        """

        res = self.results["bin"][bin_num]

        # Get the data we want to plot:
        #  for the spectral fit
        lamb = res["sl_spec"][:,0]
        obs = np.ma.masked_array(res["sl_spec"][:,1], 
                                 mask=res["sl_spec"][:,1] < 0)
        syn = res["sl_spec"][:,2]
        resid = obs - syn
        err = np.ma.masked_array(1/res["sl_spec"][:,3], 
                                 mask=res["sl_spec"][:,3] < 0)
        masked = np.ma.masked_array(resid, mask=res["sl_spec"][:,3] != 0)
        clipped = np.ma.masked_array(resid, mask=res["sl_spec"][:,3] > -1)
        #  and for mass and light fractions
        b = res["bases"]
        Z = np.sort(np.unique(b[:,5]))
        zages = []
        lightweights = []
        massweights = []
        for _Z in Z:
            _b = b[b[:,5] == _Z]
            zages.append(np.log10(_b[:,4]))
            lightweights.append(_b[:,1]) # x_j
            massweights.append(_b[:,3]) # M_cor #FIXME M_ini?????

        plt.close("all")
        # Plot the spectrum and continuum fit
        slfig = plt.figure()
        gs = gridspec.GridSpec(4, 2, height_ratios=[1, 1, 1, 1],
                               width_ratios=[3,1])
        axol = slfig.add_subplot(gs[:3,0])
        axol.plot(lamb, obs, c=OBSCOL, ls="-", lw=1, label="observed")
        axol.plot(lamb, syn, c=SYNTHCOL, ls="--", lw=1, label="starlight fit")
        axol.plot(lamb, err, c=ERRCOL, ls="-", lw=1, label="error")
        axol.set_ylabel("F$_\\lambda$ [normalised]")
        axol.set_ylim(0,2.2)
        plt.setp(axol.get_xticklabels(), visible=False) # we share the x axis

        axresid = slfig.add_subplot(gs[3,0], sharex=axol)
        axresid.plot(lamb, resid, c=OBSCOL, ls="-", lw=1, label="residual")
        axresid.plot(lamb, masked, c=MASKCOL, ls="-", lw=1, label="masked")
        axresid.plot(lamb, clipped, c=CLIPPEDCOL, ls="-", lw=1, label="clipped")
        axresid.set_xlim(np.min(lamb),np.max(lamb))
        axresid.set_xlabel("Wavelength "
                           "[{}]".format(self.lamb_units.to_string()))
        axresid.set_ylabel("Residual")
        axresid.set_yticks([-0.2,0.0,0.2,0.4,0.6])
        axresid.set_ylim(-0.3,0.7)


        # Plot the contribution of the SSPs to the population (hideously
        # written!)
        axlight = slfig.add_subplot(gs[:2,1])
        axmass = slfig.add_subplot(gs[2:,1], sharex=axlight)
        bottoml, bottomm = np.zeros(len(zages[0])), np.zeros(len(zages[0]))
        for i,_Z in enumerate(Z):
            axlight.bar(zages[i], lightweights[i], color=ZCOL[i],
                        bottom=bottoml, label="Z$_\\star$ = "+str(_Z),
                        edgecolor="none", align="center",width=0.2)
            axmass.bar(zages[i], massweights[i], color=ZCOL[i], bottom=bottomm,
                       label="Z$_\\star$ = "+str(_Z), edgecolor="none",
                       align="center",width=0.2)
            bottoml += lightweights[i]
            bottomm += massweights[i]
        axlight.plot(zages[0],[77]*len(zages[0]), marker="|", markersize=12,
                     color=SYNTHCOL, ls="none")
        lgnd = axlight.legend(loc=9, bbox_to_anchor=(0.5, 1.23),
                                frameon=False, fontsize=10, ncol=2)
        axlight.set_ylabel("x$_\\textrm{{j}}$ [\%]")
        axlight.set_ylim(0,80)
        plt.setp(axlight.get_xticklabels(), visible=False) # we share the x axis


        axmass.set_yscale("log")
        axmass.set_ylim(0.1,105)
        yformatter = ticker.FuncFormatter(
            lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),
                                                                0)))).format(y))
        axmass.yaxis.set_major_formatter(yformatter)
        axmass.set_xticks([6,7,8,9,10])
        axmass.set_xlim(5.5,10.5)
        axmass.set_xlabel("$\\log_{10}$ Age$_\\star$ [years]")
        axmass.set_ylabel("Mcor$_\\textrm{{j}}$ [\%]")

        slfig.text(0.15, 1.0, "$\chi^{{2}}/\\textrm{{dof}} = {:4.2f}$, "
                   "N$_\\textrm{{spax}} = {:3}$, A$_\\textrm{{v}} = "
                   "{:5.3f}$".format(res["chi2/Nl_eff"],len(res["x_spax"]),
                                     res["AV_min"]),color="k", size=12)
        slfig.text(0.15, 0.96, "$\\sigma_\\star = {:6.2f}$ km s$^{{-1}}$, "
                   "$v_\\star = {:6.2f}$ km s$^{{-1}}$".format(
                       res["vd_min"], res["v0_min"]), color="k", size=12)
        slfig.tight_layout()
        slfig.subplots_adjust(hspace=0.1)
        slfig.savefig(self.base_name+"_sl_fit_{}.png".format(bin_num),
                      bbox_inches="tight", dpi=300, additional_artists=(lgnd,))


    def plot_worst_fits(self, N=5):
        """
        Find the N bins with the worst chi2 from the starlight and emission
        line fits.  Plot both fits for these bins.

        """
        d = self.results["bin"]
        worst_sl_bins = sorted(d.iterkeys(), 
                               key=(lambda key:
                                    d[key]["chi2/Nl_eff"]))[-N:]
        worst_el_bins = sorted(d.iterkeys(), 
                               key=(lambda key:
                                    d[key]["emline_model"].chi2dof))[-N:]
        worst_bins = set(worst_sl_bins + worst_el_bins)
        for wb in worst_bins:
            self.plot_sl_results(wb)
            self.plot_emission_lines(wb)


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

        young = np.empty(self.data_cube.shape[1:]) * np.nan
        inter = np.empty(self.data_cube.shape[1:]) * np.nan
        old =  np.empty(self.data_cube.shape[1:]) * np.nan
        age1_str = str(int(age1/1e6))+" Myr" if age1<1e9 else \
                   str(int(age1/1e9))+" Gyr"
        age2_str = str(int(age2/1e6))+" Myr" if age2<1e9 else \
                   str(int(age2/1e9))+" Gyr"
        for b in self.results["bin"]:
            r = self.results["bin"][b]
            yngidx = r["bases"][:,4] <= age1
            intidx = (age1 < r["bases"][:,4]) & (r["bases"][:,4] < age2)
            oldidx = r["bases"][:,4] >= age2
            ynglightfrac = np.sum(r["bases"][yngidx,1])
            intlightfrac = np.sum(r["bases"][intidx,1]) 
            oldlightfrac = np.sum(r["bases"][oldidx,1])
            young[r["y_spax"],r["x_spax"]] = ynglightfrac
            inter[r["y_spax"],r["x_spax"]] = intlightfrac
            old[r["y_spax"],r["x_spax"]] = oldlightfrac

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
        fig.set_size_inches(5.,15.)

        plt.savefig(self.base_name+"_yio.pdf", bbox_inches="tight")
        print("plot saved to {}".format(self.base_name+"_yio.pdf"))    

    def plot_kinematic(self, norm_v0=0):
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
        print "kinematic_plots()"

        v0 = np.empty(self.data_cube.shape[1:]) * np.nan
        vd = np.empty(self.data_cube.shape[1:]) * np.nan
        if norm_v0 == "nucleus":
            bn = self.get_bin_num(self.nucleus)
            norm_v0 = self.results["bin"][bn]["v0_min"]
        elif not isinstance(norm_v0, (float,int)):
            print "norm_v0 must be 'nucleus' or a float/int"
            return
        for b in self.results["bin"]:
            r = self.results["bin"][b]
            v0[r["y_spax"],r["x_spax"]] = r["v0_min"] - norm_v0
            vd[r["y_spax"],r["x_spax"]] = r["vd_min"]

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

    def run_emission_lines(self, bin_num=None, vd_init=[10.,40.,70.,100.],
                           v0_init=[-300,-150,-50,0,50,150,300],
                           amp_init=[0.1,1.,10.], stddev_bounds=[5.,120.],
                           offset_bounds=[-500.,500.], weights=True):
        """
        Fit emission lines in the continuum-subtracted spectra with gaussians.

        Based on the contents of `data/emission_lines.json` (NOT AMENDABLE AT
        PRESENT!), a series of gaussians will be fit to those emission lines.
        The fitting is performed in a brute force manner over a range of 
        initial guesses to attempt to find the global minimum. As such, 
        increasing the length of the initial guesses for the width, mean and 
        amplitude (``vd_init``, ``v0_init`` and ``amp_init`` respectively) will 
        dramatically increase computing time.

        Parameters
        ----------
        bin_num : (None, int or array_like), optional
            The bin numbers to fit (defaults to None, this will
            fit all bins)
        vd_init : list
            The intial guess(es) for the velocity dispersion of the lines in 
            km/s. If None will use initial guesses of 10, 40, 70 and 100 km/s.
        v0_init : list or "vstar"
            The intial guess(es) for the velocity offset of the lines in
            km/s. The special case "vstar" will use the velocty offset
            determined by STARLIGHT and construct 6 guesses over the range
            +/-200 km/s of that value, but limited to within ``offset_bounds``
            of zero offset.
        amp_init : list
            The initial guess(es) for the amplitudes of the lines in units of 
            ``fobs_norm`` (see starlight manual).
        stddev_bounds : list
            Bounding limits for the sigma widths of the lines in km/s.
        offset_bounds : list
            Bounding limits for the offset of the mean of the lines in km/s.
        weights : bool
            Whether to include weights (1/stddev**2) in the fitting and chi2
            determination.
        """
        #TODO will need to estimate S/N of lines using noise underlying the
        #     line positions and set any S/N < 3, say, to upper limits for
        #     purposes of Z calculations etc.
        if bin_num is None:
            bin_num = self.bin_nums
        elif bin_num not in self.bin_nums:
            ValueError("bin number {} does not exist".format(bin_num))
        elif isinstance(bin_num, (int,float)):
            bin_num = [bin_num]
        print("fitting emission lines to {} bins...".format(len(bin_num)))

        # multiprocessing params for emission line fitting pool of workers
        p = mp.Pool(self.n_cpu)
        emline_results = p.map(fit_emission_lines, 
                               zip(bin_num,
                                   [self.results["bin"][bn] for bn in bin_num],
                                   repeat(self.emission_lines),
                                   repeat(vd_init),
                                   repeat(v0_init),
                                   repeat(amp_init),
                                   repeat(stddev_bounds),
                                   repeat(offset_bounds),
                                   repeat(weights)))
        p.close()
        p.join()
        print
        print("emission line fitting complete")
        for bn, emline_model in emline_results:
            # if we have no param uncertainties then the fit is at the bounds
            # of our limits and the fit is probably bad
            bad = 0
            if np.all(np.isnan(emline_model.param_uncerts)):
                bad = 1
            self.emline_output[int(bn)] = [emline_model, emline_model.chi2dof,
                                           bad]
            self.results["bin"][bn]["emline_model"] = emline_model


    def get_emission_line_fluxes(self, bin_num=None):
        """
        Calculate the equivalent widths of the emission lines.

        The values are calculated for each bin in ``bin_num`` and are based
        on the emission line models produced by run_emission_line(). Follows
        the formal definintion of emission=negative_ew, absoption=positive_ew). 

        Parameters
        ----------
        bin_num : (None, int or array_like), optional
            The bin numbers to get equivalent widths for (defaults to None, 
            this will fit all bins)
        """
        if bin_num is None:
            bin_num = self.bin_nums
        elif bin_num not in self.bin_nums:
            ValueError("Bin number {} does not exist".format(bin_num))
        elif isinstance(bin_num, (int,float)):
            bin_num = [bin_num]

        n_bins = len(bin_num)
        #FIXME current link between submodel designation and the value
        # in param_uncerts - not general!!!!
        desname = {"_0":0,"_1":3,"_2":4,"_3":5,"_4":8,"_5":9,"_6":10,"_7":11}
        for i,bn in enumerate(bin_num,1):
            print("getting emission line fluxes {}/{}\r".format(i, n_bins)),
            sys.stdout.flush()
            ew = []
            ew_uncert = []
            fluxes = []
            bin_res = self.results["bin"][bn]
            emline_model = bin_res["emline_model"]
            for submodel_name  in emline_model.submodel_names:
                #FIXME #FIXME This is implemented awfully!!!
                #FIXME these should be attributes of the model directly
                if submodel_name[-2:] in ["_0", "_1"]:
                    # i.e. Balmer
                    mean_uncert = emline_model.param_uncerts[1]
                    stddev_uncert = emline_model.param_uncerts[2]
                else:
                    # i.e. Forbidden
                    mean_uncert = emline_model.param_uncerts[6]
                    stddev_uncert = emline_model.param_uncerts[7]
                # FIXME use submodel name to determine which model uncert param
                # we need for the amplitude...
                amp_uncert = emline_model.\
                             param_uncerts[desname[submodel_name[-2:]]]
                #FIXME needs generalising!!!
                submodel = emline_model[submodel_name]
                _ew = (2*pi)**0.5 * submodel.stddev * submodel.amplitude
                # Uncertainty formula from LG IDL
                _ew_uncert = ((2*pi)**0.5 
                              * ((stddev_uncert * submodel.amplitude)**2
                                 + (submodel.stddev * amp_uncert)**2)**0.5) 
                ew.append(_ew)
                ew_uncert.append(_ew_uncert)
            bin_res["emline_ew"] = np.array(ew)
            bin_res["emline_ew_uncert"] = np.array(ew_uncert)
            bin_res["emline_fluxes"] = np.array(ew) * bin_res["fobs_norm"]
            bin_res["emline_fluxes_uncert"] = (np.array(ew_uncert) *
                                               bin_res["fobs_norm"])

    def get_emission_line_metallicities(self, bin_num=None):
        """
        Calculates several metallicity values based on popular indicators.

        Uses emission line fluxes from ``get_el_fluxes()`` and calculates PP04,
        Marino+13, Dopita+16 ... values for metallicity of each bin. These are
        saved under a new dictionary key 'metallicity' within the bins' results
        dict. All resultsare given as :math:`12+\log(O/H)`

        Parameters
        ----------
        bin_num : (None, int or array_like), optional
            The bin numbers to calculate metallicities for (defaults to None, 
            this will fit all bins)
        """

        if bin_num is None:
            bin_num = self.bin_nums
        elif bin_num not in self.bin_nums:
            ValueError("bin number {} does not exist".format(bin_num))
        elif isinstance(bin_num, (int,float)):
            bin_num = [bin_num]

        #FIXME generalise! currently line fluxes are in order of:
        # Ha, Hb, [NII]6548, [NII]6583, [SII]6716, [SII]6731, 
        # [OIII]4959, [OIII]5007
        n_bins = len(bin_num)
        for i,bn in enumerate(bin_num,1):
            print("getting emission line metallicities {}/{}\r".format(i,
                                                                       n_bins)),
            sys.stdout.flush()
            bin_res = self.results["bin"][bn]
            f = bin_res["emline_fluxes"]
            f_uncert = bin_res["emline_fluxes_uncert"]
            # We only consider fluxes with S/N > 3
            f_snr = f/f_uncert
            f[f_snr < 3] = np.nan
            fluxes = {"Halpha":f[0], "Hbeta":f[1], "[NII]6548":f[2], 
                         "[NII]6583":f[3], "[SII]6716":f[4], "[SII]6731":f[5], 
                         "[OIII]4959":f[6], "[OIII]5007":f[7]}
            bin_res["metallicity"] = {}

            #TODO allow custom indicators to be passed, or add
            # data/metallicity_indicators.json
            O3 = np.log10(fluxes["[OIII]5007"]/fluxes["Hbeta"])
            N2 = np.log10(fluxes["[NII]6583"]/fluxes["Halpha"])
            y = (np.log10(fluxes["[NII]6583"]/(fluxes["[SII]6716"] 
                                               + fluxes["[SII]6731"]))
                 + 0.264 * N2) # Dopita+16

            bin_res["metallicity"]["PP04_N2"] = 8.90 + 0.57 * (N2)
            bin_res["metallicity"]["PP04_O3N2"] = 8.73 - 0.32 * (O3/N2)
            bin_res["metallicity"]["M13"] = 8.533 - 0.214 * (O3/N2)
            bin_res["metallicity"]["D16"] = 8.77 + y

    def plot_emission_lines(self, bin_num):
        """
        Plot the residual emission line spectrum and results of emission
        line fitting for ``bin_num``

        """
        res = self.results["bin"][bin_num]

        lamb = res["sl_spec"][:,0]
        mask = res["spec"][:,3] == 2

        emline_obs = np.ma.masked_array((res["sl_spec"][:,1] -
                                         res["sl_spec"][:,2]) *
                                        res["fobs_norm"], mask=mask)
        emline_uncert = np.ma.masked_array(res["spec"][:,2], mask=mask)
        x = np.arange(lamb[0], lamb[-1], 0.25)
        emline_model = res["emline_model"](x) * res["fobs_norm"]

        plt.close("all")
        elfig,(ax1,ax2,ax3) = plt.subplots(1, 3, sharey=True,figsize=(10,5))

        for ax in (ax1,ax2,ax3):
            ax.plot(lamb, emline_obs, c=OBSCOL,lw=1)
            ax.fill_between(lamb, emline_obs-emline_uncert,
                            emline_obs+emline_uncert, color=ERRCOL)
            ax.plot(x, emline_model, c=MASKCOL,lw=1)
            ax.axvline(res["emline_model"].mean_0, ls="--", c="grey", lw=0.5)
            ax.axvline(res["emline_model"].mean_1, ls="--", c="grey", lw=0.5)
            ax.axvline(res["emline_model"].mean_2, ls="--", c="grey", lw=0.5)
            ax.axvline(res["emline_model"].mean_3, ls="--", c="grey", lw=0.5)
            ax.axvline(res["emline_model"].mean_4, ls="--", c="grey", lw=0.5)
            ax.axvline(res["emline_model"].mean_5, ls="--", c="grey", lw=0.5)
            ax.axvline(res["emline_model"].mean_6, ls="--", c="grey", lw=0.5)
            ax.axvline(res["emline_model"].mean_7, ls="--", c="grey", lw=0.5)
        ax2.set_xlabel("Wavelength [{}]".format(self.lamb_units.to_string()))
        ax1.set_ylabel("F$_\\lambda$")
        ax1.set_xlim(4830,5020)
        ax2.set_xlim(6510,6600)
        ax3.set_xlim(6690,6750)
        ax1.yaxis.tick_left()
        ax2.yaxis.set_ticks_position('none')
        ax3.yaxis.tick_right()
        ax2.xaxis.set_major_locator(ticker.MultipleLocator(20))
        ax3.xaxis.set_major_locator(ticker.MultipleLocator(20))
        ax1.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax3.spines['left'].set_visible(False)
    
        # get visible limits of data and adjust y axis limits based on those
        #FIXME relies on hardcoded emissionline limits - need to generalise
        visible = ((lamb >= 4830) & (lamb <= 5020)) | \
                  ((lamb >= 6510) & (lamb <= 6600)) | \
                  ((lamb >= 6690) & (lamb <= 6750))
        plt.ylim(np.min((emline_obs-emline_uncert)[visible]),
                 np.max((emline_obs+emline_uncert)[visible]))

        elfig.text(0.15, 1.0, "Bin num = {}, $\\bar{{x}} = {:.2f}$, "
                   "$\\bar{{y}} = {:.2f}$, $\\chi^2/\\textrm{{dof}} = {:.3f}$"
                   .format(bin_num, res["x_bar"], res["y_bar"], 
                           res["emline_model"].chi2dof),color="k", size=12)
        elfig.subplots_adjust(wspace=0.1)
        elfig.savefig(self.base_name+"_el_fit_{}.png".format(bin_num),
                      bbox_inches="tight",dpi=300)


    def plot_metallicity(self, indicator="D16", Zmin=None, Zmax=None):
        """
        Plot the metallicity of the host.

        A map of the host is shown with bins colour coded by their metallicity
        in the chosen indicator. Nan or negative values are shown as white. The
        plot is saved with suffix `_metallicity_`\ ``indicator``\ `.pdf`.

        Parameters
        ----------
        indicator: str
            The metallicity indicator to plot, current options are "PP04_N2",
            "PP04_O3N2", "M13", "D16".
        Zmin, Zmax: float, optional
            Specify the limits (in :math:`12+\log(O/H)`) for the colourmap,
            otherwise the min and max metallicities will be the bounds.
        """
        print "metallicity_plot()"

        valid_indicators = ["PP04_N2", "PP04_O3N2", "M13", "D16"]
        if indicator not in valid_indicators:
            raise AttributeError("`indicator` must be one of {}"
                                 .format(valid_indicators))

        Z = np.empty(self.data_cube.shape[1:]) * np.nan
        for b in self.results["bin"]:
            r = self.results["bin"][b]
            Z[r["y_spax"],r["x_spax"]] = r["metallicity"][indicator]

        plt.close("all")
        plt.imshow(Z, origin="lower", interpolation="none", vmin=vmin,
                   vmax=vmax)
        plt.colorbar(label="$Z$ [$12 + \log_{10}(\\textrm{O}/\\textrm{H})$]",
                     orientation="horizontal").ax.tick_params(labelsize=16)
        plt.plot(self.nucleus[0], self.nucleus[1], "kx", markersize=10)

        plt.savefig(self.base_name+"_metallicity_{}.pdf".format(indicator),
                    bbox_inches="tight")

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
                print "{} exists and clobber is false".format(outfile)
                return

        # make a dummy cube to hold the emission line spectra
        emission_line_cube = np.empty(self.data_cube.shape) * np.nan
        for b in self.results["bin"]:
            r = self.results["bin"][b]
            # subtract the synthetic spectrum from the data and account for
            # normalisation in starlight
            emission_line_spec = ((r["sl_spec"][:,1] - r["sl_spec"][:,2]) *
                                  r["fobs_norm"])
            for x, y in zip(r["x_spax"],r["y_spax"]):
                emission_line_cube[:,y,x] = emission_line_spec

        hdulist = fits.HDUList()
        hdulist.append(fits.ImageHDU(data = emission_line_cube,
                                     header=self.data_cube.header))
        hdulist.append(fits.ImageHDU(data = self.stddev_cube.data,
                                     header=self.stddev_cube.header))
        hdulist.writeto(outfile, clobber=clobber)

    def make_continuum_cube(self):
        """
        Subtract the emission line model from each bin to produce a continuum
        cube.

        To be implemented.
        """
        pass

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
                cls.__dict__ =  pickle.load(pkl_data)
            # Let's update n_cpu incase we're loading on a different machine
            cls.n_cpu = int(min(mp.cpu_count()-1,mp.cpu_count()*0.9))
            return cls
        except IOError:
            raise IOError("Couldn't create instance from pickle file"
                          " {}".format(pkl_file))
        print "loaded pkl file {}".format(pkl_file)

    def save_pkl(self, pkl_file=None):
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
        """
        if pkl_file is None:
            pkl_file = self.base_name+".pkl"
        # If the cube was compressed, we can't pickle it so need to rewrite
        # it to uncompressed format
        if isinstance(self.data_cube, fits.CompImageHDU):
            self.data_cube = fits.ImageHDU(self.data_cube.data,
                                           self.data_cube.header)
            self.stddev_cube = fits.ImageHDU(self.stddev_cube.data,
                                            self.stddev_cube.header)
        # We write to a temporary file for safety - if there's pickling
        # errors, we don't want to overwrite our previous pickle.
        temp = tempfile.mkstemp(prefix="ifuanal_", suffix=".pkl", dir=".")[1]
        print "writing to temporary pickle file {}".format(temp)
        with open(temp, "wb") as output:
            pickle.dump(self.__dict__, output)
        print "moving to {}".format(pkl_file)
        shutil.move(temp, pkl_file)
        
class EmissionLineModel(models.Gaussian1D + models.Gaussian1D +
                        models.Gaussian1D + models.Gaussian1D +
                        models.Gaussian1D + models.Gaussian1D +
                        models.Gaussian1D + models.Gaussian1D):
    # Used to explicitly define the emission line compond model in the module,
    # otherwise dill/pickle face difficulties relating to namespace issues
    pass

                
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
    vor_sn : float, optional
        The target signal-to-noise of the voronoi binning algorithm.
    sl_dir : None or str, optional
        The directory containing starlight files and bases. The default
        ``None`` will use the `starlight/` subdirectory.


    References
    ----------
    .. [SF11] Schlafly, E & Finkbeiner, D, "Measuring Reddening with Sloan\
       Digital Sky Survey Stellar Spectra and Recalibrating SFD", ApJ, 2011

    """
    def __init__(self, muse_cube, redshift, ebv="IRSA", RV=3.1, vor_sn=20,
                 sl_dir=None):
        cube_hdu = fits.open(muse_cube)
        base_name = os.path.splitext(muse_cube)[0]
        if base_name.endswith(".fits"):
            # in case original file was, e.g., *.fits.fz
            base_name = base_name[:-5] 
        data_cube = cube_hdu[1]
        stddev_cube = cube_hdu[2]
        # The MUSE STAT cube extension is the variance of the data, we want the
        # standard deviation
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            stddev_cube.data = np.sqrt(stddev_cube.data)

        if ebv == "IRSA":
            ebv = get_IRSA_ebv(data_cube.header)
        data_cube.header["IFU_EBV"] = float(ebv)

        # Add in the redshift of the target since MUSE doesn't have this in the
        # header
        data_cube.header["IFU_Z"] = float(redshift)

        super(MUSECube, self).__init__(data_cube, stddev_cube, base_name, RV,
                                       vor_sn, sl_dir)


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
    print "for location '{:}' found E(B-V) of {:5.3f}".format(coo.to_string(),
                                                              ebv)
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
    # x in CCM is in 1/microns
    x = 1/(lamb * lamb_unit.to("micron"))

    a, b = np.zeros(x.shape, x.dtype), np.ndarray(x.shape, x.dtype)

    if any((x<0.3) | (8<x)):
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
                         2.28305, 1.41338, 0), x[opt ]-1.82)

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
    Alamb = AV * (a + b/RV)

    return Alamb


def get_line_map(data_cube, lamb, line_mean, filter_width=60, cont_width=60,
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
    line_mean : float
        The central wavelength of the emission line to extract (in units of
    ``data_cube``'s header.
    filter_width : float
        The filter width to extract around the emission line.
    cont_width : float
        The width of the wavelength window to extract either side of the filter
    to define the continuum.
    cont_pad :
        The padding to apply between the edge of the filter and the start of
    continuum window.

    Returns
    -------
    2d_maps : list
        A length-4 list of 2D maps of the [blue continuum, red continuum,
        filter, line only (i.e. filter - continuum)], respectively.
    """

    # Find edges of filter and continuum windows
    low_filt = line_mean - filter_width/2.
    upp_filt = line_mean + filter_width/2.
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

    # Make maps of the three wavelength windows by averaging
    # over the spectral axis
    filter_map = np.average(data_cube.data[idx_low_filt:idx_upp_filt, :, :],
                          axis=0)
    cont1_map = np.average(data_cube.data[idx_low_cont1:idx_upp_cont1, :, :],
                          axis=0)
    cont2_map = np.average(data_cube.data[idx_low_cont2:idx_upp_cont2, :, :],
                          axis=0)

    # Determine the mean wavelength of the continuum windows
    cont1_mean = np.average((low_cont1, upp_cont1))
    cont2_mean = np.average((low_cont2, upp_cont2))
    # Interpolate their values to estimate continuum at line
    cont_interp = interp1d([cont1_mean, cont2_mean], [cont1_map, cont2_map],
                           axis=0)
    cont_at_line_mean = cont_interp(line_mean)
    # Subtract this continuum
    line_map = filter_map - cont_at_line_mean

    return [cont1_map, cont2_map, filter_map, line_map]



def resample_base(basefile, obs_lamb, delta_lamb, buffer_lamb=500):
    """
    Resample starlight base files to match ``delta_lamb`` and restrict coverage
    to limits of ``obs_lamb`` with an additional buffer.
    """
    b = np.genfromtxt(basefile)

    # Extend the ob_lamb array by buffer_lamb on each side, preserving the
    # sampling
    low = obs_lamb[0] - np.arange(1,buffer_lamb/delta_lamb)[::-1] * delta_lamb
    upp = obs_lamb[-1] + np.arange(1,buffer_lamb/delta_lamb) * delta_lamb

    interp_lamb = np.hstack((low, obs_lamb, upp))
    interp_flux = np.interp(interp_lamb, b[:,0], b[:,1])
    b_new = np.vstack((interp_lamb, interp_flux)).T
    np.savetxt(basefile, b_new)

def fit_starlight(fargs):
    """
    Perform the fitting of a spectrum with starlight and return the bin number
    and the output staright file.
    """
    bin_num, spec_file, lamb, delta_lamb, tmp_dir, bases, low_sn, upp_sn = fargs

    print "Fitting bin number {:>5}\r".format(bin_num),
    sys.stdout.flush()

    out_file = os.path.basename(spec_file)+"_out"
    # Write a grid file for starlight (see manual)
    with tempfile.NamedTemporaryFile(prefix="grid_", dir=tmp_dir, 
                                     delete=False) as grid_file:
        # we need to specifically set the precision and manually calculate the
        # wavelength window for starlight to model - in some cases the rounding
        # precision losses means our last observed wavelength is not included
        delta_lamb_syn = round(delta_lamb,8)
        lamb_syn_ini = round(lamb[0],8)
        lamb_syn_fin = lamb_syn_ini + (len(lamb) - 1) * delta_lamb_syn
        # the header
        grid_file.write("\n".join(["1", os.path.join(tmp_dir,bases,""), 
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
    sub = subprocess.Popen("nice {0} < {1} > /dev/null"\
                           .format(starlightexe,grid_file.name), shell=True,
                           cwd=tmp_dir)
    sub.wait()
    out_file_path = os.path.join(tmp_dir, out_file)
    try:
        open(out_file_path)
    except IOError:
        out_file = None
    else: 
        out_file = out_file_path
    return bin_num, out_file, spec_file


def fit_emission_lines(fargs):
    bin_num, bin_res, el, vd_init, v0_init, amp_init, stddev_b, off_b, w = fargs
    print "Fitting bin number {:>5}\r".format(bin_num),
    sys.stdout.flush()

    spec = bin_res["sl_spec"]

    if w:
        # The starlight output doesn't give the uncertainty on regions masked in
        # the fitting (i.e. emission lines!) so we need to copy these over from
        # our original spectrum and account for the normalisation of starlight
        stddev_cube = bin_res["spec"][:,2]
        stddev_cube[stddev_cube == 0] = np.nan
        spec[:,3] = bin_res["fobs_norm"]/stddev_cube**2
    else:
        spec[:,3] = np.ones(len(spec[:,0]))

    # FIXME needs generalising
    el_init = EmissionLineModel().rename("el_init_{}".format(bin_num))

    # Fix the peaks of the other gaussians relative to Halpha and [NII]
    el_init.mean_1.tied = lambda x: el["Hbeta"][0] * x.mean_0/el["Halpha"][0]
    el_init.mean_2.tied = lambda x: el["[NII]"][0] * x.mean_3/el["[NII]"][1]
    el_init.mean_4.tied = lambda x: el["[SII]"][0] * x.mean_3/el["[NII]"][1]
    el_init.mean_5.tied = lambda x: el["[SII]"][1] * x.mean_3/el["[NII]"][1]
    el_init.mean_6.tied = lambda x: el["[OIII]"][0] * x.mean_3/el["[NII]"][1]
    el_init.mean_7.tied = lambda x: el["[OIII]"][1] * x.mean_3/el["[NII]"][1]
    # Tie widths of balmer/forbidden lines separately
    el_init.stddev_1.tied  = lambda x: (x.stddev_0/x.mean_0) * x.mean_1 
    el_init.stddev_2.tied  = lambda x: (x.stddev_3/x.mean_3) * x.mean_2 
    el_init.stddev_4.tied  = lambda x: (x.stddev_3/x.mean_3) * x.mean_4 
    el_init.stddev_5.tied  = lambda x: (x.stddev_3/x.mean_3) * x.mean_5 
    el_init.stddev_6.tied  = lambda x: (x.stddev_3/x.mean_3) * x.mean_6 
    el_init.stddev_7.tied  = lambda x: (x.stddev_3/x.mean_3) * x.mean_7

    # Add some bounds that are quite generous. If the best fit is at one of
    # these limits we won't get the covariance matrix and cannot estimate
    # uncertainties. However, in this case the fit is likely to be suspect.
    el_init.mean_0.bounds = [el["Halpha"][0] * (1+off_b[0]/ckms), 
                             el["Halpha"][0] * (1+off_b[1]/ckms)]
    el_init.stddev_0.bounds =  [el["Halpha"][0] * stddev_b[0]/ckms, 
                                el["Halpha"][0] * stddev_b[1]/ckms]
    el_init.mean_3.bounds = [el["[NII]"][1] * (1+off_b[0]/ckms), 
                             el["[NII]"][1] * (1+off_b[1]/ckms)]
    el_init.stddev_3.bounds =  [el["[NII]"][1] * stddev_b[0]/ckms, 
                                el["[NII]"][1] * stddev_b[1]/ckms]

    # Mask regions away from emission lines
    lamb_mask_rows = (spec[:,0] < 4840) | \
                     ((4880 < spec[:,0]) & (spec[:,0] < 4940)) | \
                     ((5025 < spec[:,0]) & (spec[:,0] < 6530)) | \
                     ((6600 < spec[:,0]) & (spec[:,0] < 6695)) | \
                     (spec[:,0] > 6750)
    el_spec = spec[~lamb_mask_rows]

    # Use the initial guesses for param combinations later
    if v0_init == "vstar":
        offset_low = 1 + max((bin_res["v0_min"]-200), off_b[0])/ckms
        offset_high = 1 + min((bin_res["v0_min"]+200), off_b[1])/ckms
        offset_init = np.linspace(offset_low, offset_high, 6) 
    else:
        offset_init = 1 + np.array(np.array(v0_init)/ckms)
    stddev_init = np.asarray(np.array(vd_init)/ckms)

    # Make combinations of all parameter initial guesses.
    param_comb = list(product(amp_init, offset_init, stddev_init)) 
    Ncomb = len(param_comb)
    #TODO this also returns the indices of the parameters in the initial list -
    #can maybe use somehow for the fluxes/ew measurements?
    dof = len(fitting._model_to_fit_params(el_init)[0])
    chi2dof = 1e50
    best_fit = None
    # Perform minimisation with LevMar fitter for each combo to find ~global
    # minimum
    for i,comb in enumerate(param_comb, 1):
        # Make an initialisation of the model using this combo
        amp, off, sd = comb
        el_init.parameters = [amp, el["Halpha"][0] * off, el["Halpha"][0] * sd,
                              amp, el["Hbeta"][0] * off, el["Hbeta"][0] * sd,
                              amp, el["[NII]"][0] * off, el["[NII]"][0] * sd,
                              amp, el["[NII]"][1] * off, el["[NII]"][1] * sd,
                              amp, el["[SII]"][0] * off, el["[SII]"][0] * sd,
                              amp, el["[SII]"][1] * off, el["[SII]"][1] * sd,
                              amp, el["[OIII]"][0] * off, el["[OIII]"][0] * sd,
                              amp, el["[OIII]"][1] * off, el["[OIII]"][1] * sd]
        # Create the fitter and perform the fit for this combination
        levmar_fitter = fitting.LevMarLSQFitter()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            el_fit = levmar_fitter(el_init, el_spec[:,0], 
                                   el_spec[:,1]-el_spec[:,2],
                                   weights = el_spec[:,3], 
                                   maxiter=3000)
        diff = (el_spec[:,1] - el_spec[:,2]) - el_fit(el_spec[:,0])
        _chi2 = np.nansum((diff * el_spec[:,3])**2)
        _chi2dof = _chi2/dof

        if _chi2dof < chi2dof:
            chi2dof = _chi2dof
            best_fit = el_fit
            # Here we explicitly pass the fit info dict to our fitted model
            best_fit.fit_info = levmar_fitter.fit_info

    # We cannot put constraints on the parameters and still get a covariance
    # matrix, so here we correct any negative amplitudes since these are suppose
    # to be emission features (#TODO generalise this and do the opposite for
    # absorption features)
    for submodel in best_fit.submodel_names:
        # TODO add here if submodel has 'emission' in the name then min 0 (and
        # vice verse) need to name the submodels in the EmissionLineModel class
        if best_fit[submodel].amplitude.value < 0:
            best_fit[submodel].amplitude.value = 0

    best_fit.chi2dof = chi2dof
    if best_fit.fit_info["param_cov"] is None:
        print("no covariance matrix computed for bin {}, cannot"
              " compute fit uncertainties".format(bin_num))
        best_fit.param_uncerts = np.empty(dof) * np.nan
    else:
        best_fit.param_uncerts = np.diag(best_fit.fit_info["param_cov"])**0.5

    return bin_num, best_fit

def parse_starlight(sl_output):
    """
    Returns a dictionary of useful quantities from a starlight output file.
    """
    results_dict = {}
    try:
        lines = open(sl_output).readlines()
    except IOError:
        print "cannot read {}!".format(sl_output)
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
                       "vd_min", "AV_min","chi2/Nl_eff","adev",
                       "NOl_eff", "Nl_eff", "N_base","Nl_obs"]:
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
         
  
def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
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
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
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

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap     

if __name__ == "__main__":
    #TODO argparse for cli?
    pass