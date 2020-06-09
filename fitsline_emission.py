"""
Create line emission plots of FITS datacubes.

Written by Elliott Williams https://github.com/willisbillis/yt-github

Version==3.0
"""
import sys
from bisect import bisect_left
import numpy as np
from astropy.io import fits
from lmfit.models import GaussianModel
import matplotlib.pyplot as plt
plt.switch_backend('agg')

################################################################################
## CONSTANTS ##
INTERACTIVE_MODE = True # set constants here that you don't want to specify
                        # from the command line each time
# set usr_constants variables if interactive mode is off
FITSFILE = "gass_78_-48_1522043360.fits"
X, Y = 100, 30
Z_MIN_PRE = -350            # km/s - min velocity bin (before adjusting to channels)
Z_MAX_PRE = -100            # km/s - max velocity bin (before adjusting to channels)
USER_CONSTANTS = [INTERACTIVE_MODE, FITSFILE, X, Y, Z_MIN_PRE, Z_MAX_PRE]
BLOCKOUT_DATA = True        # If parts of your fits file contain no data. Set function
                            #      to determine which parts in blockout() function.
#####################
## File input/output
QUICK = 150
PRO = 1250
DPI_LEVEL = QUICK           # resolution level of output png file (dots per inch)
################################################################################
## FUNCTIONS ##
def take_closest(my_list, my_number):
    """
    Assumes myList is sorted. Returns closest value to myNumber.
    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(my_list, my_number)
    if pos == 0:
        return my_list[0]
    if pos == len(my_list):
        return my_list[-1]
    before = my_list[pos - 1]
    after = my_list[pos]
    if after - my_number < my_number - before:
        return after
    return before

def create_zlist(ffl):
    """
    Create initial velocity bins list using data from FITS data cube
    and grab FITS data cube.
    """
    # load in fits data
    hdulist = fits.open(ffl)
    hdu = hdulist[0]
    ffldata = hdu.data

    # assert that the datacube is in PPV space with units deg-deg-m/s
    hdu.header['CUNIT1'] = 'deg'
    hdu.header['CUNIT2'] = 'deg'
    hdu.header['CUNIT3'] = "m/s"

    x_dim = hdu.header['NAXIS1'] # x (deg)
    y_dim = hdu.header['NAXIS2'] # y (deg)
    z_dim = hdu.header['NAXIS3'] # z (m/s)
    print("FITS data cube dims (PPV): ", x_dim, y_dim, z_dim)

    delta_z = round(hdu.header['CDELT3'], 5)
    z_min = round(hdu.header['CRVAL3'], 5)
    z_max = round(z_min + (z_dim - 1) * delta_z, 5)

    # convert velocities from m/s to km/s
    z_list = [round(i/1000, 5) for i in np.linspace(z_min, z_max, num=z_dim)]
    return z_min, z_max, z_list, ffldata

def create_vel_list(z_min_pre, z_max_pre, z_list):
    """ Create customized velocity bins list to match channels of observations. """
    z_min = min(z_list)
    z_max = max(z_list)
    delta_z = z_list[1] - z_list[0]
    if z_min_pre < z_min or z_max_pre > z_max:
        print("Requested velocity range is outside of bounds ({0:.1f}, {1:.1f})"\
              " km/s".format(z_min, z_max))
        sys.exit()
    z_min = take_closest(z_list, z_min_pre)
    z_max = take_closest(z_list, z_max_pre)
    if z_min == z_max:
        print("Requested velocity range must be greater than {0:.4f} km/s".format(delta_z))
        sys.exit()
    if z_min > z_max:
        print("Requested velocity range is invalid ({0:.1f}, {1:.1f}) km/s. "\
              "Try switching bounds.".format(z_min, z_max))
        sys.exit()

    z_i = z_list.index(z_min)           # index of min velocity
    z_f = z_list.index(z_max)           # index of max velocity
    vel_list = np.linspace(z_min, z_max, num=z_f - z_i + 1) # list of velocity bins

    return z_i, z_f, vel_list

def inputs(usr_constants):
    """
    Used to specify variables for each run. If Interactive Mode == True, then
    variables are specified from command line prompts.
    """
    int_mode = usr_constants[0]

    if int_mode is True:
        fitsfile = input("enter a FITS filename: ").strip()
        x_coord = int(input("Input x coordinate: "))
        y_coord = int(input("Input y coordinate: "))

        z_min, z_max, z_list, ffl_data = create_zlist(fitsfile)
        vel_range = input("Current velocity range is {0:.1f} to {1:.1f} km/s. Would " \
                              "you like to change this? (y/N) ".format(z_min/1000.0, z_max/1000.0))
        if vel_range == "y":
            z_min_pre = float(input("Updated min velocity (km/s): "))
            z_max_pre = float(input("Updated max velocity (km/s): "))
        else:
            z_min_pre = round(z_min/1000.0, 5)
            z_max_pre = round(z_max/1000.0, 5)
        z_i, z_f, vel_list = create_vel_list(z_min_pre, z_max_pre, z_list)
    else:
        fitsfile = usr_constants[1]
        x_coord = usr_constants[2]
        y_coord = usr_constants[3]
        z_min, z_max, z_list, ffl_data = create_zlist(fitsfile)
        z_min_pre = usr_constants[4]
        z_max_pre = usr_constants[5]
        z_i, z_f, vel_list = create_vel_list(z_min_pre, z_max_pre, z_list)

    return ffl_data, x_coord, y_coord, z_i, z_f, vel_list

def blockout(x_coord, y_coord):
    """ Use this function to define regions of your fitsfile to avoid. """
    if (y_coord-10) > 142.0/186.0*x_coord:   # specific to GASS_***3360.fits
        return True
    return False

def line_emission(ffl_data, x_coord, y_coord, z_i, z_f, vel_list):
    """ Create line emission data and fit to gaussian from FITS data cube."""
    emission_signal = []
    for z_coord in range(z_i, z_f+1):
        emission = ffl_data[z_coord, y_coord, x_coord]
        if np.isnan(emission):
            emission_signal.append(0)
        else:
            emission_signal.append(emission)
    emission_signal = np.array(emission_signal)

    # Fit emission of LOS to gaussian curve
    mod = GaussianModel()
    pars = mod.guess(emission_signal, x=vel_list)
    out = mod.fit(emission_signal, pars, x=vel_list)
    result = np.array(out.params)
    sigma = result[0]
    center = result[1]

    print("Velocity Dispersion of peak: {:.2f} km/s".format(sigma))
    print("Velocity Center of peak: {:.2f} km/s".format(center))
    return vel_list, emission_signal, out

def plot_emission_line(vel_list, signal_list, gaussian_fit, x_coord, y_coord):
    """ Plot emission line data and gaussian best fit. """
    plt.plot(vel_list, signal_list, label="Emission line")
    plt.plot(vel_list, gaussian_fit.best_fit, label="Gaussian fit")
    plt.legend()
    plt.xlabel("Velocity (km/s)")
    plt.ylabel("K")
    plt.savefig("fits_emission_({0},{1}).png".format(x_coord, y_coord), dpi=DPI_LEVEL)
################################################################################
## PROCESSING ##

FFL_DATA, X_COORD, Y_COORD, Z_I, Z_F, VEL_LIST = inputs(USER_CONSTANTS)
VEL_LIST, EMISSION_LIST, GAUSSIAN_FIT = line_emission(FFL_DATA, X_COORD, Y_COORD,
                                                      Z_I, Z_F, VEL_LIST)
plot_emission_line(VEL_LIST, EMISSION_LIST, GAUSSIAN_FIT, X_COORD, Y_COORD)
