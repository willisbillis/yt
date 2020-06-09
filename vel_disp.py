"""
Create velocity dispersion maps of FLASH simulation data
with the post-processing software yt.

Version == 3.0

Written by Elliott Williams https://github.com/willisbillis/yt-github
"""

import time
import os
import glob
import random
import logging
import ctypes
from itertools import repeat
from math import pi, sqrt, e, log
from multiprocessing import Pool, cpu_count, Array, RLock, Process, Manager, log_to_stderr, set_start_method
from scipy.integrate import quad
from lmfit.models import GaussianModel
from mpl_toolkits.axes_grid1 import make_axes_locatable
import yt
import h5py
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')

################################################################################
## CONSTANTS ##
INTERACTIVE_MODE = True # set constants here that you don't want to specify
                        # from the command line each time
LOGGING = False          # Turn on/off logging for parallel processing
# set usr_constants variables if interactive mode is off
FULL = "n"
FL_NM = "IVC_Fall003Xv70Zv70Z10NRH43hdf5_plt_cnt_0012"
LEVEL = 1
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
        full = input("run all files? (y/N): ")
        fl_nm = input("enter a filename: ").strip()
        level = int(input("resolution level (0-4): ").strip())
        xyview = input("XY view? (y/N): ")
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

    return fixed_res_datacube, los_velocity, trans_matrix, domain, dx_resolution

def HI_filter(temp):
    """ Filter out cells that are too hot/cold
        to be picked up by HI observations """
    return HI_MIN <= temp <= HI_MAX

def model_fit(bin_weights):
    """ From an array of weights, fit a gaussian model
        and extract the std dev as the velocity dispersion """

    mod = GaussianModel()
    pars = mod.guess(bin_weights.flatten(), x=BINS[0][:-1])
    out = mod.fit(bin_weights.flatten(), pars, x=BINS[0][:-1])
    result = np.array(out.params)
    avg_vel = result[1]
    diff_sum = 0
    for vel_i, weight_i in zip(BINS[0][:-1], bin_weights.flatten()):
        diff_sum += weight_i * (vel_i-avg_vel)**2
    linefit_sigma = sqrt(diff_sum/sum(bin_weights))
    return linefit_sigma, avg_vel

def vel_disp_LOS(q_coord2, q3_dims, temp_arr, vel_arr, dens_arr, vol_arr):
    """ Create weighted array of simulated emission line along one line of sight (LOS) """
    weights = np.array([])
    weights = np.resize(weights, (NUM_BINS-1, 1))
    for q_coord3 in range(q3_dims):
        cell_temp = temp_arr[q_coord2][q_coord3]
        if HI_filter(cell_temp):
            # velocity dispersion calculation
            vel = vel_arr[q_coord2][q_coord3]
            mass = M_AVG*dens_arr[q_coord2][q_coord3]*vol_arr[q_coord2][q_coord3]
                # getting rid of thermal broadening component
            #thermal_broadening = sqrt(2 * K_BOLTZ * cell_temp / M_AVG) / 1000
            #thermal_broadening = sqrt(2 * K_BOLTZ * 100 / M_AVG) / 1000
            thermal_broadening = sqrt(2 * K_BOLTZ * 1e-1 / M_AVG) / 1000

            gaussian_mu = vel
            gaussian_sigma = thermal_broadening
            gaussian_amp = mass

            for los_cell in range(0, len(BINS[0])-1):
                emission_integral = quad(gaussian_curve, BINS[0][los_cell],
                                         BINS[0][los_cell+1],
                                         args=(gaussian_mu, gaussian_sigma, gaussian_amp))
                weights[los_cell] += emission_integral[0]

            if (gaussian_amp-sum(weights))/gaussian_amp > 0.05:
                if vel < VEL_MIN:
                    print(" Values out of range. Lower velocity limits.\n"\
                          " Out of range vel: {!s} km/s".format(vel))
                elif vel > VEL_MAX:
                    print(" Values out of range. Raise velocity limits.\n"\
                          " Out of range vel: {!s} km/s".format(vel))
    return weights

def vel_disp_row(q_coord1, shared_dict, d_map, c_map, ctr):
    """
    Calculate velocity dispersion for every line of sight (LOS) in one slab
    of the fixed_res_datacube.
    """

    num_rows = shared_dict["num_rows"]
    q2_dims = shared_dict["q2_dims"]
    q3_dims = shared_dict["q3_dims"]
    trans_matrix = shared_dict["trans_matrix"]

    temperature = np.frombuffer(shared_dict["temp"]).reshape(shared_dict["cube_dims"])
    temperature = temperature.transpose(trans_matrix)[q_coord1, :, :]
    velocity = np.frombuffer(shared_dict["velocity"]).reshape(shared_dict["cube_dims"])
    velocity = velocity.transpose(trans_matrix)[q_coord1, :, :]
    H_nuclei_density = np.frombuffer(shared_dict["Hdens"]).reshape(shared_dict["cube_dims"])
    H_nuclei_density = H_nuclei_density.transpose(trans_matrix)[q_coord1, :, :]
    cell_volume = np.frombuffer(shared_dict["cellvol"]).reshape(shared_dict["cube_dims"])
    cell_volume = cell_volume.transpose(trans_matrix)[q_coord1, :, :]

    disp_row = []
    cent_row = []

    for q_coord2 in range(q2_dims):
        weights = vel_disp_LOS(q_coord2, q3_dims, temperature, velocity,
                               H_nuclei_density, cell_volume)
        if sum(weights) != 0:
            sigma, mean_vel = model_fit(weights)
        else:
            # all cells on sightline have T outside HI_min < T < HI_max
            sigma = 0
            mean_vel = 0

        if FWHM:
            disp_row.append(2*sqrt(2*log(2))*sigma)
        else:
            disp_row.append(sigma)
        cent_row.append(mean_vel)

    ctr.value += 1
    progress = ctr.value*100.0/num_rows
    if progress % 5 < 100.0/num_rows:
        print("{0:.1f} % complete of data map".format(progress))

    d_map[q_coord1] = disp_row
    c_map[q_coord1] = cent_row

def create_data_map(file_name, level, order):
    # normal == x,z,y (0,2,1)
    # xy == y,x,z (1,0,2)
    t1 = order.index('x')
    t2 = order.index('y')
    t3 = order.index('z')
    trans_matrix = [t1,t2,t3]

    LOS_velocity = "velocity_" + order[2]
    ds = yt.load(file_name)
    all_data_level_x = ds.covering_grid(level=level,left_edge=ds.domain_left_edge,dims=ds.domain_dimensions*2**level)
    disp_array = []
    ds.periodicity = (True,True,True)
    k = 1.380649e-23 # J/K
    m_avg = 1.6737236e-27 # kg (mass of hydrogen atom)
    for q1 in range(all_data_level_x.ActiveDimensions[t1]):
        vbin = []
        for q2 in range(all_data_level_x.ActiveDimensions[t2]):
            weights = np.array([])
            weights = np.resize(weights, (NUM_BINS-1, 1))
            for q3 in range(all_data_level_x.ActiveDimensions[t3]):
                coords = tuple([q1,q2,q3][i] for i in trans_matrix)
                cell_temp = float(all_data_level_x["temperature"][coords])
                if HI_filter(cell_temp):
                    # velocity dispersion calculation
                    vel = float(all_data_level_x[LOS_velocity][coords].in_units("km/s"))
                    mass = float(m_avg*all_data_level_x["H_nuclei_density"][coords]*all_data_level_x["cell_volume"][coords])

                    # changing thermal broadening component
                    #tb = float(sqrt(2*k*cell_temp/m_avg)/1000)
                    tb = float(sqrt(2*k*100/m_avg)/1000)
                    #tb = float(sqrt(2*k*1e-1/m_avg)/1000)

                    mu = vel
                    sigma = tb
                    A = mass

                    for n in range(len(BINS[0])-1):
                        I = quad(gaussian_curve, BINS[0][n], BINS[0][n+1], args=(mu, sigma, A))
                        weights[n] += I[0]
                    if (A-sum(weights))/A > 0.05:
                        if vel < VEL_MIN:
                            print(" Values out of range. Lower velocity limits.\n"\
                                  " Out of range vel: {!s} km/s".format(vel))
                        elif vel > VEL_MAX:
                            print(" Values out of range. Raise velocity limits.\n"\
                                  " Out of range vel: {!s} km/s".format(vel))

            if sum(weights) == 0: # all cells on sightline have T outside HI_min < T < HI_max
                sigma = 0
            elif sum(weights) != 0:
                mod = GaussianModel()
                pars = mod.guess(weights.flatten(), x=BINS[0][:-1])
                out = mod.fit(weights.flatten(), pars, x=BINS[0][:-1])
                result = np.array(out.params)
                avg_vel = result[1]
                diff_sum = 0
                for a,b in zip(BINS[0][:-1], weights.flatten()):
                    diff_sum += b*(a-avg_vel)**2
                sigma = sqrt(diff_sum/sum(weights))

            if FWHM:
                vbin.append(2*sqrt(2*log(2))*sigma)
            else:
                vbin.append(sigma)
        disp_array.append(vbin)
        prog = (q1+1)*100/float(all_data_level_x.ActiveDimensions[t1])
        if prog % 5 < 100/float(all_data_level_x.ActiveDimensions[t1]):
            print("{0:.1f} % complete of data map ".format(prog))

    da = np.array(disp_array)
    left = ds.domain_left_edge.in_units("kpc")[t2]
    right = ds.domain_right_edge.in_units("kpc")[t2]
    top = ds.domain_left_edge.in_units("kpc")[t1]
    bottom = ds.domain_right_edge.in_units("kpc")[t1]
    domain = [left, right, top, bottom]
    dx = float(all_data_level_x["dx"][0, 0, 0])
    return da, domain, dx

def ang_resolution(dist, cell_width):
    """
    Calulate simulated angular resolution from width of cell in map and distance
    to target. This assumes cell width has units cm and dist has units kpc.
    """
    beam_width = (cell_width/(dist*3.08568E18))*3600
    if beam_width >= 3600:
        beam_width = beam_width/3600
        print("At distance of {} kpc, beam width is {:.1f} deg.".format(dist, beam_width))
    elif beam_width >= 60:
        beam_width = beam_width/60
        print("At distance of {} kpc, beam width is {:.1f} arcmin.".format(dist, beam_width))
    else:
        print("At distance of {} kpc, beam width is {:.1f} arcsec.".format(dist, beam_width))
    print("Velocity resolution: {0:.1f} km/s".format(BINS[1]))

def plotdata(data_array, domain_array, fl_nm, level, order, map_type):
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
    cbar.set_label("Velocity {} (km/s)".format(map_type), size=7)

    print("Creating figure...")
    fig_nm = 'velocity_{0}_{1}_{2}_lvl_{3}.png'.format(map_type[:4], ''.join(order[:2]), str(fl_nm)[-4:], level)
    plt.savefig(fig_nm, dpi=DPI_LEVEL, bbox_inches='tight')
    plt.close()
    print("File saved as: " + fig_nm)

def whole_map(fl_nm, level, xyz_order, map_type):
    da, domain, dx = create_data_map(fl_nm, level, xyz_order)
    ang_resolution(TARGET_DIST, dx)
    plotdata(da, domain, fl_nm, level, xyz_order, map_type)

def save_data(fl_list, level, order, map_type):
    """
    Save all created files to a separate directory inside the working directory.
    """
    view = ''.join(order[:2])
    directory = os.getcwd() + "/vel_{0}_{1}_lvl_{2}".format(map_type, view, level)
    if not os.path.exists(directory):
        os.makedirs(directory)
    for fl_nm in fl_list:
        working_dir = os.getcwd()
        source_file = working_dir + \
                      "/velocity_{0}_{1}_{2}_lvl_{3}.png".format(map_type, view, fl_nm[-4:], level)
        destination = working_dir + \
                      "/vel_{0}_{1}_lvl_{3}/velocity_{0}_{1}_{2}_lvl_{3}.png".format(map_type, view, fl_nm[-4:], level)
        if os.path.exists(source_file):
            print("moving {}...".format(source_file))
            os.rename(source_file, destination)
    print("Data saved to: {0}/vel_{1}_{2}_lvl_{3}".format(os.getcwd(), map_type, view, level))

def init_shared_arr(input_arr):
    """ Create empty shared 3D array for storing
        field information from FLASH datacube"""
    input_arr = np.array(input_arr)
    shared_array_base = Array(ctypes.c_double, input_arr.shape[0] * input_arr.shape[1]
                              * input_arr.shape[2])
    shared_array = np.ctypeslib.as_array(shared_array_base.get_obj())
    shared_array_np = np.frombuffer(shared_array).reshape(input_arr.shape)
    np.copyto(shared_array_np, input_arr)
    return shared_array

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

def map_solve(filename, level):
    ytcube = yt.load(filename)
    ytcube.periodicity = (True, True, True)
    fixed_res_datacube = ytcube.covering_grid(level=level, left_edge=ytcube.domain_left_edge,
                                              dims=ytcube.domain_dimensions*2**level)
    flat_temp = np.array(fixed_res_datacube["temperature"].flatten())
    print("Determining optimal configuration...")
    relevant_voxels = sum(map(HI_filter, flat_temp))
    print(relevant_voxels, "relevant voxels")
    if relevant_voxels <= 1e4:
        implementation = "whole_map"
    if 1e4 < relevant_voxels <= 1e5:
        implementation = "row_split"
    if relevant_voxels > 1e5:
        implementation = "LOS_split"
    return implementation
################################################################################
## PROCESSING ##
if __name__ == "__main__":
    set_start_method("spawn")
    FILE_LIST, LEVEL, XYZ_ORDER = inputs(USER_CONSTANTS)
    IMPLEMENTATION = map_solve(FILE_LIST[-1], LEVEL)
    print("Implementation:", IMPLEMENTATION)

    if LOGGING:
        LOGGER = log_to_stderr()
        LOGGER.setLevel(logging.DEBUG)

    START_TIME = time.time()

    if IMPLEMENTATION == "whole_map":
        ARG_LIST = zip(FILE_LIST, repeat(LEVEL), repeat(XYZ_ORDER), repeat("disp"))
        for arg in ARG_LIST:
            p = Process(target=whole_map, args=(arg[0], arg[1], arg[2], arg[3]))
            p.start()
            p.join()
        save_data(FILE_LIST, LEVEL, XYZ_ORDER, "disp")

    elif IMPLEMENTATION == "row_split":
        for fl in FILE_LIST:
            # Make necessary values and datasets available to all of the worker processes
            # so that each dataset doesn't need to be loaded in every time a calculation
            # is done.
            manager = Manager()
            SHARED_DICT = {}

            # Create fixed resolution data cube
            DATACUBE, LOS_VEL, T_MATRIX, DOMAIN, CELL_RES = create_fixed_res_array(fl, LEVEL, XYZ_ORDER)

            NUM_ROWS = int(DATACUBE.ActiveDimensions[T_MATRIX[0]])

            COMPLETED_ROWS = manager.Value('i', 0)
            SHARED_DISP_MAP = manager.list(np.zeros(NUM_ROWS))
            SHARED_CENT_MAP = manager.list(np.zeros(NUM_ROWS))
            SHARED_DICT["velocity"] = init_shared_arr(DATACUBE[LOS_VEL].in_units('km/s'))
            SHARED_DICT["cellvol"] = init_shared_arr(DATACUBE['cell_volume'])
            SHARED_DICT["Hdens"] = init_shared_arr(DATACUBE['H_nuclei_density'])
            SHARED_DICT["temp"] = init_shared_arr(DATACUBE['temperature'])
            SHARED_DICT["cube_dims"] = DATACUBE.ActiveDimensions
            SHARED_DICT["num_rows"] = DATACUBE.ActiveDimensions[T_MATRIX[0]]
            SHARED_DICT["q2_dims"] = DATACUBE.ActiveDimensions[T_MATRIX[1]]
            SHARED_DICT["q3_dims"] = DATACUBE.ActiveDimensions[T_MATRIX[2]]
            SHARED_DICT["trans_matrix"] = T_MATRIX

            # Create argument list to be iterated through by workers
            arg_list = zip(range(NUM_ROWS), repeat(SHARED_DICT), repeat(SHARED_DISP_MAP),
                           repeat(SHARED_CENT_MAP), repeat(COMPLETED_ROWS))

            # Send jobs out to workers
            CPU_CORES = num_cores(NUM_ROWS)
            LOCK = RLock()
            with Pool(processes=(CPU_CORES)) as POOL:
                RESULT = POOL.starmap(vel_disp_row, arg_list, chunksize=CHUNKSIZE)

            # Sort output dispersion and centers data into maps and save to h5 files
            final_disp_map = np.array(SHARED_DISP_MAP)
            final_cent_map = np.array(SHARED_CENT_MAP)
            print("Dispersion map created, saving data...")
            if os.path.isfile('disp_data_'+fl+'.h5'):
                h5py.File('disp_data_'+fl+'.h5', 'w').close() # clear file if exists
            with h5py.File('disp_data_'+fl+'.h5', 'w') as h5_store:
                h5_store.create_dataset('disp_map', data=final_disp_map)
                h5_store.create_dataset('cent_map', data=final_cent_map)

            # Create plot of velocity disperion map
            plotdata(final_disp_map, DOMAIN, fl, LEVEL, XYZ_ORDER, "dispersion")
            plotdata(final_cent_map, DOMAIN, fl, LEVEL, XYZ_ORDER, "center")
            ang_resolution(TARGET_DIST, CELL_RES)

        save_data(FILE_LIST, LEVEL, XYZ_ORDER, "disp")
        save_data(FILE_LIST, LEVEL, XYZ_ORDER, "cent")

    elif IMPLEMENTATION == "LOS_split":
        for fl in FILE_LIST:
            # Make necessary values and datasets available to all of the worker processes
            # so that each dataset doesn't need to be loaded in every time a calculation
            # is done.
            manager = Manager()
            SHARED_DICT = {}

            # Create fixed resolution data cube
            DATACUBE, LOS_VEL, T_MATRIX, DOMAIN, CELL_RES = create_fixed_res_array(fl, LEVEL, XYZ_ORDER)

            NUM_ROWS = int(DATACUBE.ActiveDimensions[T_MATRIX[0]])

            COMPLETED_ROWS = manager.Value('i', 0)
            SHARED_DISP_MAP = manager.list(np.zeros(NUM_ROWS))
            SHARED_CENT_MAP = manager.list(np.zeros(NUM_ROWS))
            SHARED_DICT["velocity"] = init_shared_arr(DATACUBE[LOS_VEL].in_units('km/s'))
            SHARED_DICT["cellvol"] = init_shared_arr(DATACUBE['cell_volume'])
            SHARED_DICT["Hdens"] = init_shared_arr(DATACUBE['H_nuclei_density'])
            SHARED_DICT["temp"] = init_shared_arr(DATACUBE['temperature'])
            SHARED_DICT["cube_dims"] = DATACUBE.ActiveDimensions
            SHARED_DICT["num_rows"] = DATACUBE.ActiveDimensions[T_MATRIX[0]]
            SHARED_DICT["q2_dims"] = DATACUBE.ActiveDimensions[T_MATRIX[1]]
            SHARED_DICT["q3_dims"] = DATACUBE.ActiveDimensions[T_MATRIX[2]]
            SHARED_DICT["trans_matrix"] = T_MATRIX

            # Create argument list to be iterated through by workers
            arg_list = zip(range(NUM_ROWS), repeat(SHARED_DICT), repeat(SHARED_DISP_MAP),
                           repeat(SHARED_CENT_MAP), repeat(COMPLETED_ROWS))

            # Send jobs out to workers
            CPU_CORES = num_cores(NUM_ROWS)
            LOCK = RLock()
            with Pool(processes=(CPU_CORES)) as POOL:
                RESULT = POOL.starmap(vel_disp_row, arg_list, chunksize=CHUNKSIZE)

            # Sort output dispersion and centers data into maps and save to h5 files
            final_disp_map = np.array(SHARED_DISP_MAP)
            final_cent_map = np.array(SHARED_CENT_MAP)
            print("Dispersion map created, saving data...")
            if os.path.isfile('disp_data_'+fl+'.h5'):
                h5py.File('disp_data_'+fl+'.h5', 'w').close() # clear file if exists
            with h5py.File('disp_data_'+fl+'.h5', 'w') as h5_store:
                h5_store.create_dataset('disp_map', data=final_disp_map)
                h5_store.create_dataset('cent_map', data=final_cent_map)

            # Create plot of velocity disperion map
            plotdata(final_disp_map, DOMAIN, fl, LEVEL, XYZ_ORDER, "dispersion")
            plotdata(final_cent_map, DOMAIN, fl, LEVEL, XYZ_ORDER, "center")
            ang_resolution(TARGET_DIST, CELL_RES)

        save_data(FILE_LIST, LEVEL, XYZ_ORDER, "disp")
        save_data(FILE_LIST, LEVEL, XYZ_ORDER, "cent")
    STOP_TIME = time.time()
    print(timer(STOP_TIME-START_TIME))
