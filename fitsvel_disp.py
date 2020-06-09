"""
Create velocity dispersion and centers maps of FITS datacubes.

Written by Elliott Williams https://github.com/willisbillis/yt-github

Version==3.0
"""
import os
import time
import sys
from itertools import repeat
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
INTERACTIVE_MODE = True     # set constants here that you don't want to specify
                            # from the command line each time
## File input/output
FITSFILE = "gass_78_-48_1522043360.fits"
DISP_OUTPUT_FL = "fits_vel_disp.png"
CENT_OUTPUT_FL = "fits_disp_cents.png"
QUICK = 150
PRO = 1250
DPI_LEVEL = QUICK           # resolution level of output png file (dots per inch)

## Velocity dispersion parameters
SIG_THRESHOLD = 0.053       # set to rms noise level of detections. used to
                            #      filter out non-signals in gaussian fitting
Y_INIT = 0
Z_MIN_PRE = -350            # km/s - min velocity bin (before adjusting to channels)
Z_MAX_PRE = -100            # km/s - max velocity bin (before adjusting to channels)
DISP_MAX = 30               # km/s - largest expected value for dispersion, used to
                            #      filter out poorly fit emission lines
BLOCKOUT_DATA = True        # If parts of your FITS file contain no data. Set function
                            #      to determine which parts in blockout() function.
USER_CONSTANTS = [INTERACTIVE_MODE, FITSFILE, Z_MIN_PRE, Z_MAX_PRE]
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
    coord_sys = WCS(hdu.header, naxis=2)

    x_dim = hdu.header['NAXIS1'] # x (deg)
    y_dim = hdu.header['NAXIS2'] # y (deg)
    z_dim = hdu.header['NAXIS3'] # z (m/s)
    print("FITS data cube dims (PPV): ", x_dim, y_dim, z_dim)

    delta_z = round(hdu.header['CDELT3'], 5)
    z_min = round(hdu.header['CRVAL3'], 5)
    z_max = round(z_min + (z_dim - 1) * delta_z, 5)

    # convert velocities from m/s to km/s
    z_list = [round(i/1000, 5) for i in np.linspace(z_min, z_max, num=z_dim)]
    return z_min, z_max, z_list, coord_sys, ffldata

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
    Used to specify variables for each run. If Interactive Mode is True, then
    variables are specified from command line prompts.
    """
    int_mode = usr_constants[0]

    if int_mode is True:
        fitsfile = input("enter a FITS filename: ").strip()

        z_min, z_max, z_list, coord_sys, ffl_data = create_zlist(fitsfile)
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
        z_min, z_max, z_list, coord_sys, ffl_data = create_zlist(fitsfile)
        z_min_pre = usr_constants[2]
        z_max_pre = usr_constants[3]
        z_i, z_f, vel_list = create_vel_list(z_min_pre, z_max_pre, z_list)

    x_dim = ffl_data.shape[2]
    y_dim = ffl_data.shape[1]
    return x_dim, y_dim, z_i, z_f, vel_list, coord_sys, ffl_data

def blockout(x_coord, y_coord):
    """ Use this function to define regions of your fitsfile to avoid. """
    if (y_coord-10) > 142.0/186.0*x_coord:   # specific to GASS_***3360.fits
        return True
    return False

def fit_line_emission(x_coord, y_coord, z_i, z_f, vel_bins, datacube):
    """ Do a gaussian fit to a single emission line to find its center
    and dispersion (sigma)"""
    emission_signal = []
    for z_coord in range(z_i, z_f+1):
        emission = datacube[z_coord, y_coord, x_coord]
        if np.isnan(emission):
            emission_signal.append(0.0)
        else:
            emission_signal.append(emission)
    emission_signal = np.array(emission_signal)

    # Fit emission of LOS to gaussian curve
    mod = GaussianModel()
    pars = mod.guess(emission_signal, x=vel_bins)
    out = mod.fit(emission_signal, pars, x=vel_bins)
    result = np.array(out.params)
    sigma = result[0]
    center = result[1]
    height = result[4]
    return sigma, center, height

def fit_row(x_dims, y_coord, z_i, z_f, vel_bins, datacube):
    """ Fits emission lines for an entire row of coordinates. Because this is a
    computationally expensive process, the image is divided into chunks of rows. """

    print("Running y = {0}/{1}".format(y_coord+1, datacube.shape[1]))
    vel_max = max(vel_bins)
    vel_min = min(vel_bins)
    x_vals = list(range(x_dims))
    sigma_list = list(np.zeros(x_dims))
    center_list = list(np.zeros(x_dims))
    no_signal = []
    case2 = []

    for x_coord in x_vals:
        if BLOCKOUT_DATA and blockout(x_coord, y_coord):
            pass
        else:
            sigma, center, height = fit_line_emission(x_coord, y_coord, z_i, z_f,
                                                      vel_bins, datacube)
            # Sort results of gaussian fitting into positive/negative signals
            if height > SIG_THRESHOLD and vel_min < center < vel_max:
                if sigma < DISP_MAX:
                    sigma_list[x_coord] = sigma
                    center_list[x_coord] = center
                else:
                    case2.append((x_coord, y_coord))
            elif height < SIG_THRESHOLD or center < vel_min or center > vel_max:
                no_signal.append((x_coord, y_coord))

    h5_sigma_list = [y_coord] + sigma_list      # add y index value to store in h5 data file
    np_center_list = [y_coord] + [center_list]  # same as above but for 2D numpy array
    return h5_sigma_list, no_signal, case2, np_center_list

def num_cores(num_jobs):
    """ Determine the number of workers needed in multiprocessing Pool """
    try:
        cpu_cores = int(os.environ["PBS_NP"])
    except KeyError:
        print("Cluster network not detected. Using multiprocessing to detect cores.")
        cpu_cores = cpu_count() - 1
    if cpu_cores > num_jobs:
        cpu_cores = num_jobs
    return cpu_cores
################################################################################
## PROCESSING ##
if __name__ == '__main__':
    START_TIME = time.time()

    # Initialize necessary variables for computation
    X_DIM, Y_DIM, Z_I, Z_F, VEL, MAP_WCS, DATA = inputs(USER_CONSTANTS)


    # Run parallel processes
    CPU_CORES = num_cores(Y_DIM)
    ARG_INPUTS = zip(repeat(X_DIM), range(Y_DIM), repeat(Z_I), repeat(Z_F), repeat(VEL), repeat(DATA))
    with Pool(processes=CPU_CORES) as POOL:
        DATA_MASTER = POOL.starmap(fit_row, ARG_INPUTS)

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
    print("time1:", timer(STOP_TIME1-START_TIME))

    # Use data from hdf5 file to make plot of dispersion map
    if os.path.isfile('disp_data.h5'):
        with h5py.File('disp_data.h5', 'r') as hf:
            DATA = hf['alldata'][:]
            DATA = DATA[DATA[:, 0].argsort()] #reorder data by first entry of each element
            DATA = np.array([list(i[1:]) for i in DATA]) # take off ordering and
                                                         # return np.array
            FIG = plt.figure()
            AXES = WCSAxes(FIG, [0.1, 0.1, 0.8, 0.8], wcs=MAP_WCS)
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
    print("time2:", timer(STOP_TIME2-START_TIME))

    FIG = plt.figure()
    AXES = WCSAxes(FIG, [0.1, 0.1, 0.8, 0.8], wcs=MAP_WCS)
    FIG.add_axes(AXES)
    IMAGE = plt.imshow(CENTER_MAIN, origin="lower", cmap="jet", vmin=min(VEL), vmax=max(VEL))

    plt.xlabel("longitude")
    plt.ylabel("latitude")
    CBAR = plt.colorbar(pad=0)
    CBAR.ax.tick_params(labelsize=8)
    CBAR.set_label("Dispersion Centers (km/s)", size=7)
    plt.savefig(CENT_OUTPUT_FL, dpi=DPI_LEVEL, bbox_inches='tight')
    plt.cla()
