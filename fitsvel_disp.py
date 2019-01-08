import os
import h5py
import time
import sys
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
from bisect import bisect_left
from lmfit.models import GaussianModel
from multiprocessing import Pool, cpu_count
from setup import timer

################################################################################
## CONSTANTS ##
Interactive_Mode = False
# set constants here that you don't want to specify from the command line each time

## File input/output
fitsfile = "gass_78_-48_1522043360.fits"
disp_output_fl = "fits_vel_disp.png"
cent_output_fl = "fits_disp_cents.png"
fig_quality = 750

## Parallel processing parameters
chunksize = 3 # best if a multiple of processes for efficiency
maxtasks = chunksize*5 # max tasks per child before printing results and starting new child process

## Dispersion parameters
sig_threshold = 0.053 # set to rms noise level of detections. used to filter out non-signals in gaussian fitting
y_init = 0
z_min_pre = -350 # km/s - min velocity bin before adjusting to channels
z_max_pre = -100 # km/s - max velocity bin before adjusting to channels
constants = [sig_threshold,z_min_pre,z_max_pre] # list of user defined constants to export to functions
################################################################################
## FUNCTIONS ##
def flatlist(list):
    return [item for sublist in list for item in sublist]

def takeClosest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return after
    else:
       return before

def inputs(ffl,usr_constants):
    """
    Used to specify variables for each run. If Interactive_Mode == True, then

    variables are specified from command line prompts.
    """
    hdulist = fits.open(ffl)
    hdu = hdulist[0]
    ffldata = hdu.data

    hdu.header['CUNIT1'] = 'deg'
    hdu.header['CUNIT2'] = 'deg'
    hdu.header['CUNIT3'] = "m/s"
    coord_sys = WCS(hdu.header,naxis=2)

    x_dim = hdu.header['NAXIS1']            # x (deg), 187 pixels
    y_dim = hdu.header['NAXIS2']            # y (deg), 251 pixels
    z_dim = hdu.header['NAXIS3']            # z (vel), 1201 frames
    dz = hdu.header['CDELT3']               # width of each velocity channel
    z_min = hdu.header['CRVAL3']            # minimum velocity of domain
    z_max = z_min + z_dim*dz                # maximum velocity of domain
    z_list = [i for i in np.arange(z_min,z_max,dz)] # list of velocity bins for entire domain

    if y_dim % chunksize != 0:
        y_dim = y_dim - (y_dim % chunksize)
        ffldata = ffldata[:,:y_dim,:]
        print "Y dimension not divisible by chunksize. Adjusted to {0} from original {1}.".format(y_dim,hdu.header['NAXIS2'])

    if Interactive_Mode == True:
        vel_range = raw_input("Current velocity range is {0:.1f} - {1:.1f} km/s. "
                              "Would you like to change this? (y/N) ".format(z_min/1000.0,z_max/1000.0))
        if vel_range == "y":
            z_min_pre = float(raw_input("Updated min velocity (km/s): ")) * 1000
            z_max_pre = float(raw_input("Updated max velocity (km/s): ")) * 1000
    else:
        z_min_pre = usr_constants[1] * 1000
        z_max_pre = usr_constants[2] * 1000

    if z_min_pre < z_min or z_max_pre > z_max:
        print "Requested velocity range is outside of bounds ({0:.1f},{1:.1f}) km/s".format(z_min/1000.0,z_max/1000.0)
        sys.exit()
    z_min = takeClosest(z_list,z_min_pre)
    z_max = takeClosest(z_list,z_max_pre)
    if z_min_pre != z_min or z_max_pre != z_max:
        print "Velocity bounds adjusted to ({0:.1f},{1:.1f}) km/s to match channel width.".format(z_min/1000.0,z_max/1000.0)

    z_i = z_list.index(z_min)           # index of min velocity
    z_f = z_list.index(z_max)           # index of max velocity
    vel_list = np.array([i/1000 for i in np.arange(z_min,z_max,dz)]) # list of velocity bins for chosen range

    return x_dim,y_dim,z_i,z_f,vel_list,coord_sys,ffldata

def line_emission((y,x_dims,z_i,z_f,vel,datacube)):
    print "Running y = {0}/{1}".format(y+1,datacube.shape[1])
    vel_max = max(vel)
    vel_min = min(vel)
    x_vals = [i for i in xrange(0,x_dims)]
    vbin = [0 for i in x_vals]

    no_signal = []
    case2 = []
    center_list = [0 for i in x_vals]

    for x in x_vals:
        if (y-10) > 142.0/186.0*x:               # specific to GASS_***3360.fits
            continue
        else:
            K = []
            for z in xrange(z_i,z_f):
                k = datacube[z,y,x]
                if np.isnan(k) == True:
                    K.append(0)
                else:
                    K.append(k)
            K = np.array(K)

            # Fit emission of LOS to gaussian curve
            mod = GaussianModel()
            pars = mod.guess(K, x=vel)
            out  = mod.fit(K, pars, x=vel)
            result = np.array(out.params)
            sigma = result[0]
            center = result[1]
            height = result[4]
            center_list[x] = center

            # Sort results of gaussian fitting into positive/negative signals
            if height > sig_threshold and vel_min < center < vel_max:
                if sigma < 30.0:
                    vbin[x] = sigma
                    continue
                else:
                    case2.append((x,y))
                    continue
            elif height < sig_threshold or center < vel_min or center > vel_max:
                no_signal.append((x,y))
                continue

    h5_vbin = [y] + vbin    # add y index value to vbin to store in h5 data file
    np_center_list = [y] + [center_list] # same as above but for 2D numpy array
    return h5_vbin,no_signal,case2,np_center_list
################################################################################
## PROCESSING ##

if __name__ == '__main__':
    time0 = time.time()
    cpus = cpu_count()
    pool = Pool(processes=cpus-1, maxtasksperchild=maxtasks)

    # Initialize necessary variables for computation
    x_dim, y_dim, z_i, z_f, vel, CAR_WCS, data = inputs(fitsfile, constants)

    # Run parallel processes
    arg_inputs = [(y, x_dim, z_i, z_f, vel, data) for y in range(y_init, y_dim)]
    data_master = pool.map(line_emission, arg_inputs, chunksize=chunksize)

    # Sort output data
    disp_data = np.array([i[0] for i in data_master])
    no_signal_main = flatlist([i[1] for i in data_master])
    case2_main = flatlist([i[2] for i in data_master])
    center_main = sorted([i[3] for i in data_master]) # sort by y index
    center_main = np.array([i[1] for i in center_main]) #take off ordering and return array

    # Save dispersion map to h5 file
    if os.path.isfile('disp_data.h5'):
        h5py.File('disp_data.h5', 'w').close()
    with h5py.File("disp_data.h5", 'a') as hf:
        hf.create_dataset('alldata', data=disp_data, maxshape=disp_data.shape)
    pool.close()
    pool.join()

    # Save alternate cases to npy files
    np.save("no_signal", no_signal_main)
    np.save("case2", case2_main)

    time1 = time.time()
    print "time1:", timer(time1-time0)

    # Use data from hdf5 file to make plot of dispersion map
    if os.path.isfile('disp_data.h5'):
        with h5py.File('disp_data.h5', 'r') as hf:
            data = hf['alldata'][:]
            data = data[data[:,0].argsort()] #reorder data by first entry of each element
            data = np.array([list(i[1:]) for i in data]) # take off ordering and
                                                         # return np.array
            fig = plt.figure()
            ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=CAR_WCS)
            fig.add_axes(ax)
            im = plt.imshow(data, origin = "lower", cmap='jet')

            plt.xlabel("longitude")
            plt.ylabel("latitude")
            cbar = plt.colorbar(pad = 0)
            cbar.ax.tick_params(labelsize=8)
            cbar.set_label("Velocity Dispersion (km/s)", size=7)
            plt.savefig(disp_output_fl, dpi=fig_quality, bbox_inches='tight')
            plt.cla()
    time2 = time.time()
    print "time2:", timer(time2-time0)

    fig = plt.figure()
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=CAR_WCS)
    fig.add_axes(ax)
    im = plt.imshow(center_main, origin="lower", cmap="jet", vmin=min(vel), vmax=max(vel))

    plt.xlabel("longitude")
    plt.ylabel("latitude")
    cbar = plt.colorbar(pad = 0)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label("Dispersion Centers (km/s)", size=7)
    plt.savefig(cent_output_fl, dpi=fig_quality, bbox_inches='tight')
    plt.cla()
