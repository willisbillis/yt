"""
Create velocity dispersion and centers maps of FITS datacubes.

Written by Elliott Williams https://github.com/willisbillis/yt-github

Version==2.0
"""
import os
import time
import sys
from bisect import bisect_left
from multiprocessing import Pool, cpu_count
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
from lmfit.models import GaussianModel
import h5py
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')

################################################################################
## CONSTANTS ##
INTERACTIVE_MODE = False     # set constants here that you don't want to specify
                            # from the command line each time

## File input/output
FITSFILE = "gass_78_-48_1522043360.fits"
DISP_OUTPUT_FL = "fits_vel_disp.png"
CENT_OUTPUT_FL = "fits_disp_cents.png"
QUICK = 150
PRO = 1250
DPI_LEVEL = QUICK           # resolution level of output png file (dots per inch)

## Parallel processing parameters
CHUNKSIZE = 3               # best if a multiple of processes for efficiency
MAXTASKS = CHUNKSIZE * 5    # max tasks per child before printing results and
                            #      starting new child process

## Dispersion parameters
SIG_THRESHOLD = 0.053       # set to rms noise level of detections. used to
                            #      filter out non-signals in gaussian fitting
Y_INIT = 0
Z_MIN_PRE = -350            # km/s - min velocity bin (before adjusting to channels)
Z_MAX_PRE = -100            # km/s - max velocity bin (before adjusting to channels)
DISP_MAX = 30               # km/s - largest expected value for disperion, used to
                            #      filter out poorly fit emission lines
BLOCKOUT_DATA = True        # If parts of your fits file contain no data. Set function
                            #      to determine which parts in blockout() function.
USER_CONSTANTS = [INTERACTIVE_MODE, Z_MIN_PRE, Z_MAX_PRE]
                 # list of user defined constants to export to functions
################################################################################
## FUNCTIONS ##
def flatlist(input_list):
    """ Remove inner lists within lists """
    return [item for sublist in input_list for item in sublist]

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

def timer(t_sec):
    """
    Output time in more readable units.
    """
    if t_sec <= 60:
        readout = "Elapsed time: {:.1f} sec".format(t_sec)
    elif t_sec <= 7200:
        t_min = t_sec / 60.0
        readout = "Elapsed time: {:.1f} min".format(t_min)
    return readout

def inputs(ffl, usr_constants):
    """
    Used to specify variables for each run. If Interactive Mode is True, then
    variables are specified from command line prompts.
    """
    hdulist = fits.open(ffl)
    hdu = hdulist[0]
    ffldata = hdu.data
    header = hdu.header

    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'
    header['CUNIT3'] = "m/s"
    coord_sys = WCS(header, naxis=2)

    x_dim = header['NAXIS1']                # x (deg), 187 pixels
    y_dim = header['NAXIS2']                # y (deg), 251 pixels
    z_dim = header['NAXIS3']                # z (vel), 1201 frames
    delta_z = header['CDELT3']                   # width of each velocity channel
    z_min = header['CRVAL3']                # minimum velocity of domain
    z_max = z_min + z_dim * delta_z                    # maximum velocity of domain
    z_list = list(np.arange(z_min, z_max, delta_z))  # list of velocity bins for entire domain

    if y_dim % CHUNKSIZE != 0:
        y_dim = y_dim - (y_dim % CHUNKSIZE)
        ffldata = ffldata[:, :y_dim, :]
        print "Y dimension not divisible by chunksize of {0}. Adjusted to {1}" \
              " from original {2}.".format(CHUNKSIZE, y_dim, header['NAXIS2'])

    int_mode = usr_constants[0]

    if int_mode is True:
        vel_range = raw_input("Current velocity range is {0:.1f} to {1:.1f} km/s. Would " \
                              "you like to change this? (y/N) ".format(z_min/1000.0, z_max/1000.0))
        if vel_range == "y":
            z_min_pre = float(raw_input("Updated min velocity (km/s): ")) * 1000
            z_max_pre = float(raw_input("Updated max velocity (km/s): ")) * 1000
    else:
        z_min_pre = usr_constants[1] * 1000
        z_max_pre = usr_constants[2] * 1000

    if z_min_pre < z_min or z_max_pre > z_max:
        print "Requested velocity range is outside of bounds ({0:.1f}, {1:.1f})"\
              " km/s".format(z_min/1000.0, z_max/1000.0)
        sys.exit()
    z_min = take_closest(z_list, z_min_pre)
    z_max = take_closest(z_list, z_max_pre)
    if z_min_pre != z_min or z_max_pre != z_max:
        print "Velocity bounds adjusted to ({0:.1f}, {1:.1f}) km/s to match" \
              " channel width.".format(z_min/1000.0, z_max/1000.0)

    z_i = z_list.index(z_min)           # index of min velocity
    z_f = z_list.index(z_max)           # index of max velocity
    vel_list = np.arange(z_min, z_max, delta_z) / 1000 # list of velocity bins

    return x_dim, y_dim, z_i, z_f, vel_list, coord_sys, ffldata

def blockout(x_coord, y_coord):
    """ Use this function to define regions of your fitsfile to avoid. """
    if (y_coord-10) > 142.0/186.0*x_coord:   # specific to GASS_***3360.fits
        return True
    return False

def fit_line_emission(x_coord, y_coord, z_i, z_f, vel_bins, datacube):
    """ Do a gaussian fit to a single emission line to find its center
    and dispersion (sigma)"""
    emission_line = []
    for z_coord in xrange(z_i, z_f):
        k = datacube[z_coord, y_coord, x_coord]
        if np.isnan(k):
            emission_line.append(0.0)
        else:
            emission_line.append(k)
    emission_line = np.array(emission_line)
    # Fit emission of LOS to gaussian curve
    mod = GaussianModel()
    pars = mod.guess(emission_line, x=vel_bins)
    out = mod.fit(emission_line, pars, x=vel_bins)
    result = np.array(out.params)
    sigma = result[0]
    center = result[1]
    height = result[4]
    return sigma, center, height

def fit_row((x_dims, y_coord, z_i, z_f, vel_bins, datacube)):
    """ Fits emission lines for an entire row of coordinates. Because this is a
    computationally expensive process, the image is divided into chunks of rows. """

    print "Running y = {0}/{1}".format(y_coord+1, datacube.shape[1])
    vel_max = max(vel_bins)
    vel_min = min(vel_bins)
    x_vals = list(xrange(0, x_dims))
    sigma_list = center_list = list(np.zeros(x_dims))
    no_signal = []
    case2 = []

    for x_coord in x_vals:
        if BLOCKOUT_DATA and blockout(x_coord, y_coord):
            continue
        else:
            sigma, center, height = fit_line_emission(x_coord, y_coord, z_i, z_f,
                                                      vel_bins, datacube)
            # Sort results of gaussian fitting into positive/negative signals
            if height > SIG_THRESHOLD and vel_min < center < vel_max:
                if sigma < DISP_MAX:
                    sigma_list[x_coord] = sigma
                    center_list[x_coord] = center
                    continue
                else:
                    case2.append((x_coord, y_coord))
                    continue
            elif height < SIG_THRESHOLD or center < vel_min or center > vel_max:
                no_signal.append((x_coord, y_coord))
                continue

    h5_sigma_list = [y_coord] + sigma_list      # add y index value to store in h5 data file
    np_center_list = [y_coord] + [center_list]  # same as above but for 2D numpy array
    return h5_sigma_list, no_signal, case2, np_center_list
################################################################################
## PROCESSING ##

if __name__ == '__main__':
    START_TIME = time.time()
    CPU_CORES = cpu_count()
    POOL = Pool(processes=(CPU_CORES-1), maxtasksperchild=MAXTASKS)

    # Initialize necessary variables for computation
    X_DIM, Y_DIM, Z_I, Z_F, VEL, CAR_WCS, DATA = inputs(FITSFILE, USER_CONSTANTS)

    # Run parallel processes
    ARG_INPUTS = [(X_DIM, Y, Z_I, Z_F, VEL, DATA) for Y in range(Y_INIT, Y_DIM)]
    DATA_MASTER = POOL.map(fit_row, ARG_INPUTS, chunksize=CHUNKSIZE)
    POOL.close()
    POOL.join()

    # Sort output data
    DISP_DATA = np.array([i[0] for i in DATA_MASTER])
    NO_SIGNAL_MAIN = flatlist([i[1] for i in DATA_MASTER])
    CASE2_MAIN = flatlist([i[2] for i in DATA_MASTER])
    CENTER_MAIN = sorted([i[3] for i in DATA_MASTER]) # sort by y index
    CENTER_MAIN = np.array([i[1] for i in CENTER_MAIN]) #take off ordering and return array

    # Save dispersion map to h5 file
    if os.path.isfile('disp_data.h5'):
        h5py.File('disp_data.h5', 'w').close()
    with h5py.File("disp_data.h5", 'a') as hf:
        hf.create_dataset('alldata', data=DISP_DATA, maxshape=DISP_DATA.shape)

    # Save alternate cases to npy files
    np.save("no_signal", NO_SIGNAL_MAIN)
    np.save("case2", CASE2_MAIN)

    STOP_TIME1 = time.time()
    print "time1:", timer(STOP_TIME1-START_TIME)

    # Use data from hdf5 file to make plot of dispersion map
    if os.path.isfile('disp_data.h5'):
        with h5py.File('disp_data.h5', 'r') as hf:
            DATA = hf['alldata'][:]
            DATA = DATA[DATA[:, 0].argsort()] #reorder data by first entry of each element
            DATA = np.array([list(i[1:]) for i in DATA]) # take off ordering and
                                                         # return np.array
            FIG = plt.figure()
            AXES = WCSAxes(FIG, [0.1, 0.1, 0.8, 0.8], wcs=CAR_WCS)
            FIG.add_axes(AXES)
            IMAGE = plt.imshow(DATA, origin="lower", cmap='jet')

            plt.xlabel("longitude")
            plt.ylabel("latitude")
            CBAR = plt.colorbar(pad=0)
            CBAR.ax.tick_params(labelsize=8)
            CBAR.set_label("Velocity Dispersion (km/s)", size=7)
            plt.savefig(DISP_OUTPUT_FL, dpi=DPI_LEVEL, bbox_inches='tight')
            plt.cla()
    STOP_TIME2 = time.time()
    print "time2:", timer(STOP_TIME2-START_TIME)

    FIG = plt.figure()
    AXES = WCSAxes(FIG, [0.1, 0.1, 0.8, 0.8], wcs=CAR_WCS)
    FIG.add_axes(AXES)
    IMAGE = plt.imshow(CENTER_MAIN, origin="lower", cmap="jet", vmin=min(VEL), vmax=max(VEL))

    plt.xlabel("longitude")
    plt.ylabel("latitude")
    CBAR = plt.colorbar(pad=0)
    CBAR.ax.tick_params(labelsize=8)
    CBAR.set_label("Dispersion Centers (km/s)", size=7)
    plt.savefig(CENT_OUTPUT_FL, dpi=DPI_LEVEL, bbox_inches='tight')
    plt.cla()
