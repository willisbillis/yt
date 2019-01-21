"""
Create velocity dispersion maps of FLASH simulation data

with the post-processing software yt.

Written by Elliott Williams https://github.com/willisbillis/yt-github
"""

import time
import sys
import os
import glob
import random
from math import pi, sqrt, e, log
from multiprocessing import Pool, cpu_count
from scipy.integrate import quad
from lmfit.models import GaussianModel
from mpl_toolkits.axes_grid1 import make_axes_locatable
import yt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

################################################################################
## CONSTANTS ##
INTERACTIVE_MODE = True # set constants here that you don't want to specify
                        # from the command line each time
# set usr_constants variables if interactive mode is off
FULL = "y"
FL_NM = ""
LEVEL = 0
XYVIEW = "n"
USER_CONSTANTS = [INTERACTIVE_MODE, FULL, FL_NM, LEVEL, XYVIEW]
#############

## File input/output
QUICK = 150
PRO = 1250
DPI_LEVEL = QUICK        # resolution level of output png file (dots per inch)

## Parallel processing parameters
CHUNKSIZE = 3            # best if a multiple of processes for efficiency
MAXTASKS = CHUNKSIZE * 5 # max tasks per child before printing results and
                         #      starting new child process

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
        fl_nm = os.getcwd() +  "/" + raw_input("enter a filename: ").strip()
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
        file_list = glob.glob(fl_nm[0:-4] + "*")
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

def fixed_res_array(file_name, level, order):
    """
    Create a 3D grid with equal sized cells from the input file and go through
    each line of sight (LOS) and calculate velocity dispersion. Puts dispersion
    values in 2D numpy map.
    """
    # normal == x,z,y (0,2,1)
    # xy == y,x,z (1,0,2)
    t_coord1 = order.index('x')
    t_coord2 = order.index('y')
    t_coord3 = order.index('z')
    trans_matrix = [t_coord1, t_coord2, t_coord3]

    los_velocity = "velocity_" + order[2]
    ytcube = yt.load(file_name)
    all_data_level_x = ytcube.covering_grid(level=level, left_edge=ytcube.domain_left_edge,
                                            dims=ytcube.domain_dimensions*2**level)
    disp_array = []
    ytcube.periodicity = (True, True, True)
    k = 1.380649e-23 # J/K
    m_avg = 1.6737236e-27 # kg (mass of hydrogen atom)
    for q_coord1 in xrange(all_data_level_x.ActiveDimensions[t_coord1]):
        vbin = []
        for q_coord2 in xrange(all_data_level_x.ActiveDimensions[t_coord2]):
            weights = np.array([])
            weights = np.resize(weights, (NUM_BINS-1, 1))
            for q_coord3 in xrange(all_data_level_x.ActiveDimensions[t_coord3]):
                coords = tuple([q_coord1, q_coord2, q_coord3][i] for i in trans_matrix)
                temp = float(all_data_level_x["temperature"][coords])
                if temp <= HI_MAX and temp >= HI_MIN:
                    # velocity dispersion calculation
                    cell_weights = np.array([])
                    cell_weights = np.resize(weights, (NUM_BINS-1, 1))
                    vel = float(all_data_level_x[los_velocity][coords].in_units("km/s"))
                    mass = float(m_avg * all_data_level_x["H_nuclei_density"][coords]
                                 * all_data_level_x["cell_volume"][coords])
                    thermal_broadening = float(sqrt(2 * k * temp / m_avg) / 1000)

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
                vbin.append(2*sqrt(2*log(2))*linefit_sigma)
            elif FWHM is False:
                vbin.append(linefit_sigma)
        disp_array.append(vbin)
        prog = (q_coord1+1)*100/float(all_data_level_x.ActiveDimensions[t_coord1])
        if prog % 5 < 100/float(all_data_level_x.ActiveDimensions[t_coord1]):
            print "{0:.1f} %".format(prog)

    da_npgrid = np.array(disp_array)
    print "fixed resolution array created"
    left = ytcube.domain_left_edge.in_units("kpc")[t_coord2]
    right = ytcube.domain_right_edge.in_units("kpc")[t_coord2]
    top = ytcube.domain_left_edge.in_units("kpc")[t_coord1]
    bottom = ytcube.domain_right_edge.in_units("kpc")[t_coord1]
    domain = [left, right, top, bottom]
    dx_resolution = all_data_level_x["dx"]
    return da_npgrid, domain, dx_resolution

def ang_resolution(dist, cell_width):
    """
    Calulate simulated angular resolution from width of cell in map and distance
    to target.
    """
    beam_width = (cell_width[0, 0, 0]/(dist*3.08568E18))*3600
    if beam_width >= 3600:
        beam_width = beam_width/3600
        print "At distance of {} kpc, beam width is {.1f} deg.".format(dist, beam_width)
    elif beam_width >= 60:
        beam_width = beam_width/60
        print "At distance of {} kpc, beam width is {.1f} arcmin.".format(dist, beam_width)
    else:
        print "At distance of {} kpc, beam width is {.1f} arcsec.".format(dist, beam_width)
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

    print "plot created. Saving figure..."
    fig_nm = 'velocity_disp_{0}_lvl_{1}.png'.format(str(fl_nm)[-4:], level)
    plt.savefig(fig_nm, dpi=DPI_LEVEL, bbox_inches='tight')
    plt.close()
    print "File saved as: " + fig_nm
    return

def save_data(fl_list, level, order):
    """
    Save all created files to a separate directory inside the working directory.
    """
    directory = os.getcwd() + "/vel_disp_lvl_{0}".format(level)
    if not os.path.exists(directory):
        os.makedirs(directory)
    for fl_nm in fl_list:
        working_dir = os.getcwd()
        source_file = working_dir + "/velocity_{0}disp_{1}_lvl_{2}.png".format(''.join(order[:1]),
                                                                               fl_nm[-4:], level)
        destination = working_dir + "/vel_{0}disp_lvl_{2}/" + \
                      "velocity_{0}disp_{1}_lvl_{2}.png".format(''.join(order[:1]),
                                                                fl_nm[-4:], level)
        if os.path.exists(source_file):
            print "moving {}...".format(source_file)
            os.rename(source_file, destination)
    print "Data saved to: {1}/vel_{0}disp_lvl_{2}".format(''.join(order[:1]), os.getcwd(), level)

def main((fl_nm, level, xyz_order)):
    """
    Main execution function. Create map from data cube, calculate angular
    resolution, and plot data.
    """
    data_map, domain, dx_res = fixed_res_array(fl_nm, level, xyz_order)
    ang_resolution(TARGET_DIST, dx_res)
    plotdata(data_map, domain, fl_nm, level, xyz_order)
    return
################################################################################
## PROCESSING ##
if __name__ == "__main__":
    FILE_LIST, LEVEL, XYZ_ORDER = inputs(USER_CONSTANTS)
    ARG_LIST = zip(FILE_LIST, [LEVEL]*len(FILE_LIST), [XYZ_ORDER]*len(FILE_LIST))

    START_TIME = time.time()
    CPU_CORES = cpu_count()
    POOL = Pool(processes=(CPU_CORES-1), maxtasksperchild=MAXTASKS)
    RESULT = POOL.map(main, ARG_LIST, chunksize=CHUNKSIZE)
    POOL.close()
    POOL.join()
    save_data(FILE_LIST, LEVEL, XYZ_ORDER)
    STOP_TIME = time.time()
    print timer(STOP_TIME-START_TIME)
