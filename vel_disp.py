"""
Create velocity dispersion maps of FLASH simulation data
with the post-processing software yt.

Version == 2.1

Written by Elliott Williams https://github.com/willisbillis/yt-github
"""

import time
import sys
import os
import glob
import random
import logging
from math import pi, sqrt, e, log
from multiprocessing import Pool, cpu_count, Manager, Array, Lock, log_to_stderr
import ctypes
from scipy.integrate import quad
from lmfit.models import GaussianModel
from mpl_toolkits.axes_grid1 import make_axes_locatable
import yt
import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
plt.switch_backend('agg')

################################################################################
## CONSTANTS ##
INTERACTIVE_MODE = False # set constants here that you don't want to specify
                        # from the command line each time
LOGGING = True          # Turn on/off logging for parallel processing
# set usr_constants variables if interactive mode is off
FULL = "n"
FL_NM = "dual_4r100_hdf5_plt_cnt_0000"
LEVEL = 2
XYVIEW = "n"
USER_CONSTANTS = [INTERACTIVE_MODE, FULL, FL_NM, LEVEL, XYVIEW]
#############

## File input/output
QUICK = 150
PRO = 1250
DPI_LEVEL = QUICK        # resolution level of output png file (dots per inch)

## Parallel processing parameters
CHUNKSIZE = 1            # num of jobs to submit to worker at a time

## Velocity dispersion parameters
TARGET_DIST = 65         # (kpc) Assumed distance to the target in map - used to
                         #      calculate angular resolution
HI_MIN = 0               # (K) Ionization cutoff temp for H - only material in
                         #      this temperature range will be "detected"
HI_MAX = 10000           # (K)
FWHM = False             # Turn on/off FWHM measurements (assumes vel disp along
                         #      LOS is Gaussian, which is valid)
NUM_BINS = 200           # Bins in each emission line. More bins >>
                         #      higher velocity/frequency resolution
VEL_MIN = -100           # (km/s)
VEL_MAX = 100            # (km/s)
BINS = np.linspace(VEL_MIN, VEL_MAX, num=NUM_BINS, retstep=True)

## Physics constants
K_BOLTZ = 1.380649e-23   # J/K Boltzmann constant
M_AVG = 1.6737236e-27    # kg (mass of hydrogen atom)
################################################################################
## FUNCTIONS ##
def inputs(usr_constants):
    """
    Used to specify variables for each run. If Interactive Mode == True, then
    variables are specified from command line prompts.
    """
    int_mode = usr_constants[0]

    if int_mode is True:
        full = raw_input("run all files? (y/N): ")
        fl_nm = raw_input("enter a filename: ").strip()
        level = int(raw_input("resolution level (0-4): ").strip())
        xyview = raw_input("XY view? (y/N): ")
    else:
        full = usr_constants[1]

        fl_nm = usr_constants[2]
        level = usr_constants[3]
        xyview = usr_constants[4]

    if xyview == "y":
        xyz_order = ['x', 'y', 'z']
    else:
        xyz_order = ['x', 'z', 'y']
    if full == "y":
        file_list = glob.glob(fl_nm[0:-4] + "*[0-9]")
    else:
        file_list = [fl_nm]

    return file_list, level, xyz_order

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

def gaussian_curve(x_val, mean_mu, sigma, amp):
    """
    Gaussian curve formula for fitting emission line emission.
    """
    return (amp/(sigma*sqrt(2*pi)))*e**(-0.5*((x_val-mean_mu)/sigma)**2)

def create_fixed_res_array(file_name, level, order):
    """
    Create the 3D grid with equal sized cells from the input yt file, and
    get measurements/resolution for the 3D grid domain. Also, arrange the
    transformation matrix and LOS velocity so that the velocity dispersion
    calculation is done in the correct dimension for the orientation requested.
    """
    # normal (xz) view == x,z,y (0,2,1)
    # xy view == y,x,z (1,0,2)
    t_coord1 = order.index('x')
    t_coord2 = order.index('y')
    t_coord3 = order.index('z')
    trans_matrix = [t_coord1, t_coord2, t_coord3]
    los_velocity = "velocity_" + order[2]

    ytcube = yt.load(file_name)
    ytcube.periodicity = (True, True, True)
    fixed_res_datacube = ytcube.covering_grid(level=level, left_edge=ytcube.domain_left_edge,
                                              dims=ytcube.domain_dimensions*2**level)

    left = ytcube.domain_left_edge.in_units("kpc")[t_coord2]
    right = ytcube.domain_right_edge.in_units("kpc")[t_coord2]
    top = ytcube.domain_left_edge.in_units("kpc")[t_coord1]
    bottom = ytcube.domain_right_edge.in_units("kpc")[t_coord1]
    domain = [left, right, top, bottom]
    dx_resolution = float(fixed_res_datacube["dx"][0, 0, 0])

    return fixed_res_datacube, trans_matrix, los_velocity, domain, dx_resolution

def vel_disp_row(q_coord1):
    """
    Calculate velocity dispersion for every line of sight (LOS) in one slab
    of the fixed_res_datacube.
    """
    global NUM_BINS, M_AVG, HI_MAX, HI_MIN, K_BOLTZ, VEL_MIN, VEL_MAX, BINS, FWHM

    global shared_dict

    temperature = np.frombuffer(shared_dict["temp"]).reshape(shared_dict["cube_dims"])
    velocity = np.frombuffer(shared_dict["velocity"]).reshape(shared_dict["cube_dims"])
    H_nuclei_density = np.frombuffer(shared_dict["Hdens"]).reshape(shared_dict["cube_dims"])
    cell_volume = np.frombuffer(shared_dict["cellvol"]).reshape(shared_dict["cube_dims"])

    q2_dims = shared_dict["q2_dims"]
    q3_dims = shared_dict["q3_dims"]
    trans_matrix = shared_dict["trans_matrix"]

    disp_row = []
    for q_coord2 in xrange(q2_dims):
        weights = np.array([])
        weights = np.resize(weights, (NUM_BINS-1, 1))
        for q_coord3 in xrange(q3_dims):
            idx_coord1 = [q_coord1, q_coord2, q_coord3][trans_matrix[0]]
            idx_coord2 = [q_coord1, q_coord2, q_coord3][trans_matrix[1]]
            idx_coord3 = [q_coord1, q_coord2, q_coord3][trans_matrix[2]]
            temp = float(temperature[idx_coord1][idx_coord2][idx_coord3])
            if temp <= HI_MAX and temp >= HI_MIN:
                # velocity dispersion calculation
                cell_weights = np.array([])
                cell_weights = np.resize(weights, (NUM_BINS-1, 1))
                vel = float(velocity[idx_coord1][idx_coord2][idx_coord3])
                mass = float(M_AVG * H_nuclei_density[idx_coord1][idx_coord2][idx_coord3]
                             * cell_volume[idx_coord1][idx_coord2][idx_coord3])
                thermal_broadening = float(sqrt(2 * K_BOLTZ * temp / M_AVG) / 1000)
                gaussian_mu = vel
                gaussian_sigma = thermal_broadening
                gaussian_amp = mass

                for los_cell in xrange(0, len(BINS[0])-1):
                    emission_integral = quad(gaussian_curve, BINS[0][los_cell],
                                             BINS[0][los_cell+1],
                                             args=(gaussian_mu, gaussian_sigma, gaussian_amp))
                    cell_weights[los_cell] += emission_integral[0]
                    weights[los_cell] += emission_integral[0]
                if (gaussian_amp-sum(weights))/gaussian_amp > 0.05:
                    if vel < VEL_MIN:
                        print " Values out of range. Lower velocity limits.\n"\
                              " Out of range vel: {!s} km/s".format(vel)
                    if vel > VEL_MAX:
                        print " Values out of range. Raise velocity limits.\n"\
                              " Out of range vel: {!s} km/s".format(vel)
                    sys.exit()
        if sum(weights) == 0: # if so, all cells on sightline have T outside HI_min < T < HI_max
            linefit_sigma = 0
        elif sum(weights) != 0:
            mod = GaussianModel()
            pars = mod.guess(weights.flatten(), x=BINS[0][:-1])
            out = mod.fit(weights.flatten(), pars, x=BINS[0][:-1])
            result = np.array(out.params)
            avg_vel = result[1]
            diff_sum = 0
            for vel_i, weight_i in zip(BINS[0][:-1], weights.flatten()):
                diff_sum += weight_i * (vel_i-avg_vel)**2
            linefit_sigma = sqrt(diff_sum/sum(weights))

        if FWHM is True:
            disp_row.append(2*sqrt(2*log(2))*linefit_sigma)
        elif FWHM is False:
            disp_row.append(linefit_sigma)

    shared_data_map[q_coord1] = disp_row
    return

def ang_resolution(dist, cell_width):
    """
    Calulate simulated angular resolution from width of cell in map and distance
    to target. This assumes cell width has units cm and dist has units kpc.
    """
    beam_width = (cell_width/(dist*3.08568E18))*3600
    if beam_width >= 3600:
        beam_width = beam_width/3600
        print "At distance of {} kpc, beam width is {:.1f} deg.".format(dist, beam_width)
    elif beam_width >= 60:
        beam_width = beam_width/60
        print "At distance of {} kpc, beam width is {:.1f} arcmin.".format(dist, beam_width)
    else:
        print "At distance of {} kpc, beam width is {:.1f} arcsec.".format(dist, beam_width)
    print "Velocity resolution: {0:.1f} km/s".format(BINS[1])

def plotdata(data_array, domain_array, fl_nm, level, order):
    """
    Create png plot of dispersion map with 2D numpy data array.
    """
    matplotlib.rcParams['font.size'] = 10
    fig_num = random.randint(0, 500)
    plt.figure(fig_num)
    axes = plt.gca()
    image = axes.imshow(data_array, origin="lower", aspect="equal",
                        extent=domain_array, cmap="jet", interpolation='none')
    plt.xlabel(order[1] + " (kpc)")
    plt.ylabel(order[0] + " (kpc)")
    divider = make_axes_locatable(axes)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(image, cax=cax)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label("Velocity Dispersion (km/s)", size=7)

    print "Creating figure..."
    fig_nm = 'velocity_disp{0}_{1}_lvl_{2}.png'.format(''.join(order[:2]), str(fl_nm)[-4:], level)
    plt.savefig(fig_nm, dpi=DPI_LEVEL, bbox_inches='tight')
    plt.close()
    print "File saved as: " + fig_nm
    return

def save_data(fl_list, level, order):
    """
    Save all created files to a separate directory inside the working directory.
    """
    view = ''.join(order[:2])
    directory = os.getcwd() + "/vel_disp_{0}_lvl_{1}".format(view, level)
    if not os.path.exists(directory):
        os.makedirs(directory)
    for fl_nm in fl_list:
        working_dir = os.getcwd()
        source_file = working_dir + \
                      "/velocity_disp{0}_{1}_lvl_{2}.png".format(view, fl_nm[-4:], level)
        destination = working_dir + \
                      "/vel_disp_{0}_lvl_{2}/velocity_disp{0}_{1}_lvl_{2}.png".format(view, fl_nm[-4:], level)
        if os.path.exists(source_file):
            print "moving {}...".format(source_file)
            os.rename(source_file, destination)
    print "Data saved to: {0}/vel_disp_{1}_lvl_{2}".format(os.getcwd(), view, level)

def init_shared_arr(input_arr):
    input_arr = np.array(input_arr)
    shared_array_base = Array(ctypes.c_double, input_arr.shape[0]*input_arr.shape[1]*input_arr.shape[2])
    shared_array = np.ctypeslib.as_array(shared_array_base.get_obj())
    shared_array_np = np.frombuffer(shared_array).reshape(input_arr.shape)
    np.copyto(shared_array_np, input_arr)
    return shared_array

def init_worker(lock, datacube, los_vel, t_matrix):
    lock.acquire()
    shared_dict["velocity"] = init_shared_arr(datacube[los_vel].in_units('km/s'))
    shared_dict["cellvol"] = init_shared_arr(datacube['cell_volume'])
    shared_dict["Hdens"] = init_shared_arr(datacube['H_nuclei_density'])
    shared_dict["temp"] = init_shared_arr(datacube['temperature'])
    shared_dict["cube_dims"] = datacube.ActiveDimensions
    shared_dict["q2_dims"] = datacube.ActiveDimensions[t_matrix[1]]
    shared_dict["q3_dims"] = datacube.ActiveDimensions[t_matrix[2]]
    shared_dict["trans_matrix"] = t_matrix
    lock.release()

################################################################################
## PROCESSING ##
if __name__ == "__main__":
    FILE_LIST, LEVEL, XYZ_ORDER = inputs(USER_CONSTANTS)
    try:
        CPU_CORES = int(os.environ["PBS_NP"])
    except KeyError:
        CPU_CORES = cpu_count() - 1

    if LOGGING:
        LOGGER = log_to_stderr()
        LOGGER.setLevel(logging.INFO)

    START_TIME = time.time()
    for fl in FILE_LIST:
        # Create fixed resolution data cube
        DATACUBE, T_MATRIX, LOS_VEL, DOMAIN, CELL_RES = create_fixed_res_array(fl, LEVEL, XYZ_ORDER)

        # Make necessary values and datasets available to all of the worker processes
        # so that each dataset doesn't need to be loaded in every time a calculation
        # is done.
        manager = Manager()

        COMPLETED_ROWS = manager.Value('i', 0)
        shared_dict = {}

        # Create argument list to be iterated through by workers
        NUM_ROWS = int(DATACUBE.ActiveDimensions[T_MATRIX[0]])
        arg_list = xrange(NUM_ROWS)

        # Send jobs out to workers
        shared_data_map = manager.list(np.zeros(NUM_ROWS))
        LOCK = Lock()
        POOL = Pool(processes=(CPU_CORES), initializer=init_worker,
                    initargs=(LOCK, DATACUBE, LOS_VEL, T_MATRIX))
        RESULT = POOL.map_async(vel_disp_row, arg_list, chunksize=CHUNKSIZE)
        POOL.close()
        POOL.join()
        PROGRESS = 0
        while True:
            if RESULT.ready():
                break
            REMAINING = RESULT._number_left
            PROGRESS_NEW = (NUM_ROWS - REMAINING)*100.0/NUM_ROWS
            if PROGRESS_NEW % 5 < 100.0/NUM_ROWS and PROGRESS_NEW != PROGRESS:
                print "{0:.1f} % complete of data map".format(PROGRESS_NEW)
                PROGESS = PROGRESS_NEW
            time.sleep(0.5)

        # Sort output dispersion data into map and save to h5 file
        shared_data_map = np.array(shared_data_map)
        print "Dispersion map created, saving data..."
        if os.path.isfile('disp_data_'+fl+'.h5'):
            h5py.File('disp_data_'+fl+'.h5', 'w').close() # clear file if exists
        with h5py.File('disp_data_'+fl+'.h5') as h5_store:
            h5_store.create_dataset('disp_map', data=shared_data_map)

        # Create plot of velocity disperion map
        plotdata(shared_data_map, DOMAIN, fl, LEVEL, XYZ_ORDER)
        ang_resolution(TARGET_DIST, CELL_RES)


    save_data(FILE_LIST, LEVEL, XYZ_ORDER)
    STOP_TIME = time.time()
    print timer(STOP_TIME-START_TIME)
