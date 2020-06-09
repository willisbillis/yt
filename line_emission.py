"""
Create line emission plots of FLASH datacubes.

Written by Elliott Williams https://github.com/willisbillis/yt-github

Version==3.0
"""
from math import pi, sqrt, e, log
from lmfit.models import GaussianModel
import yt
from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')

################################################################################
## CONSTANTS ##
INTERACTIVE_MODE = True # set constants here that you don't want to specify
                        # from the command line each time
# set usr_constants variables if interactive mode is off
FL_NM = "sing_hdf5_plt_cnt_0050"
LEVEL = 4
XYVIEW = "n"
X, Y = 128, 128
USER_CONSTANTS = [INTERACTIVE_MODE, FL_NM, LEVEL, XYVIEW, X, Y]
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
#####################
## File input/output
QUICK = 150
PRO = 1250
DPI_LEVEL = QUICK        # resolution level of output png file (dots per inch)
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
        fl_nm = input("enter a filename: ").strip()
        level = int(input("resolution level (0-4): ").strip())
        xyview = input("XY view? (y/N): ")
        x_coord = int(input("Input x coordinate: "))
        y_coord = int(input("Input y coordinate: "))
    else:
        fl_nm = usr_constants[1]
        level = usr_constants[2]
        xyview = usr_constants[3]
        x_coord = usr_constants[4]
        y_coord = usr_constants[5]

    if xyview == "y":
        xyz_order = ['x', 'y', 'z']
    else:
        xyz_order = ['x', 'z', 'y']

    return fl_nm, level, xyz_order, x_coord, y_coord

def gaussian_curve(x_val, mean_mu, sigma, amp):
    """ Gaussian curve formula for fitting emission line emission. """
    return (amp/(sigma*sqrt(2*pi)))*e**(-0.5*((x_val-mean_mu)/sigma)**2)

def create_fixed_res_array(file_name, level, order, x_coord, y_coord):
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
    print(fixed_res_datacube["dx"].transpose(trans_matrix).shape)

    temperature = fixed_res_datacube['temperature'].transpose(trans_matrix)[x_coord, y_coord, :]
    temperature = np.array(temperature)
    vel = fixed_res_datacube[los_velocity].transpose(trans_matrix)[x_coord, y_coord, :]
    vel = np.array(vel.in_units("km/s"))
    dens = fixed_res_datacube['H_nuclei_density'].transpose(trans_matrix)[x_coord, y_coord, :]
    dens = np.array(dens)
    vol = fixed_res_datacube['cell_volume'].transpose(trans_matrix)[x_coord, y_coord, :]
    vol = np.array(vol)

    fields_datacube = [temperature, vel, dens, vol]

    dx_resolution = float(fixed_res_datacube["dx"][0, 0, 0])

    return fields_datacube, trans_matrix, dx_resolution

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
    return linefit_sigma, avg_vel, out

def vel_disp_LOS(data_arr):
    """ Create weighted array of simulated emission line along one line of sight (LOS) """
    temp_arr, vel_arr, dens_arr, vol_arr = data_arr
    weights = np.array([])
    weights = np.resize(weights, (NUM_BINS-1, 1))
    for q_coord3 in range(len(temp_arr)):
        cell_temp = temp_arr[q_coord3]
        if HI_filter(cell_temp):
            # velocity dispersion calculation
            vel = vel_arr[q_coord3]
            mass = M_AVG*dens_arr[q_coord3]*vol_arr[q_coord3]

            # changing thermal broadening component
            thermal_broadening = sqrt(2 * K_BOLTZ * cell_temp / M_AVG) / 1000
            #thermal_broadening = sqrt(2 * K_BOLTZ * 100 / M_AVG) / 1000
            #thermal_broadening = sqrt(2 * K_BOLTZ * 1e-1 / M_AVG) / 1000

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
"""
def line_emission():

    if sum(weights) == 0: # all cells on sightline have T outside HI_min < T < HI_max
        sigma = 0
    elif sum(weights) != 0:
        mod = GaussianModel()
        pars = mod.guess(weights.flatten(), x=BINS[0][:-1])
        print(type(weights.flatten()))
        out  = mod.fit(weights.flatten(), pars, x=BINS[0][:-1])
        #print(out.fit_report(min_correl=0.25))
        result = np.array(out.params)
        mod_sigma = result[0]
        print(mod_sigma)

        avg_vel = result[1]
        diff_sum = 0
        for a, b in zip(BINS[0][:-1], weights.flatten()):
            diff_sum += b*(a-avg_vel)**2
        sigma = sqrt(diff_sum/sum(weights))
        print(sigma)

    percent_diff = 200*abs(sigma-mod_sigma)/(sigma+mod_sigma)
    print(str(percent_diff) + " %")

    print("Velocity resolution: {0:.1f} km/s".format(BINS[1]))
    print("cells measured along LOS: " + str(cl_cnt) + "/" + str(len(velz)))
"""
def plot_emission_line(vel_list, signal_list, gaussian_fit, x_coord, y_coord):
    """ Plot emission line data and gaussian best fit. """
    plt.plot(vel_list, signal_list, label="Emission line")
    plt.plot(vel_list, gaussian_fit.best_fit, label="Gaussian fit")
    plt.legend()
    plt.xlabel("Velocity (km/s)")
    plt.ylabel("Intensity")
    plt.savefig("emission_({0},{1}).png".format(x_coord, y_coord), dpi=DPI_LEVEL)

################################################################################
## PROCESSING ##
FL_NM, LEVEL, XYZ_ORDER, X_COORD, Y_COORD = inputs(USER_CONSTANTS)
DATACUBE, T_MATRIX, CELL_RES = create_fixed_res_array(FL_NM, LEVEL, XYZ_ORDER, X_COORD, Y_COORD)
WEIGHTS = vel_disp_LOS(DATACUBE)
if sum(WEIGHTS) == 0: # all cells on sightline have T outside HI_min < T < HI_max
    SIGMA = CENTER = 0
elif sum(WEIGHTS) != 0:
    SIGMA, CENTER, MODEL_FIT = model_fit(WEIGHTS)
    plot_emission_line(BINS[0][:-1], WEIGHTS, MODEL_FIT, X_COORD, Y_COORD)
if FWHM:
    SIGMA = 2*sqrt(2*log(2))*SIGMA
print("Velocity Dispersion of peak: {:.2f} km/s".format(SIGMA))
print("Velocity Center of peak: {:.2f} km/s".format(CENTER))
ang_resolution(TARGET_DIST, CELL_RES)
