"""
Create velocity dispersion maps of FLASH simulation data

with the post-processing software yt.

Written by Elliott Williams https://github.com/willisbillis/yt-github

Version==1.0
"""

import yt,time,matplotlib,sys,os,glob,random
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import multiprocessing as mp
from math import pi,sqrt,e
from scipy.integrate import quad
from lmfit.models import GaussianModel

def timer(t):
    # outputs time in more readable units
    if t <= 60:
        return "Elapsed time: %.1f sec" % t
    elif t <= 7200:
        t = float(t)/60
        return "Elapsed time: %.1f min" % t
################################################################################
# set resolution level of output png file (dots per inch)
quick = 150
pro = 1250
dpi_level = quick
# Turn on/off FWHM measurements (assumes vel disp along LOS is Gaussian, which is valid)
FWHM = False
# Ionization cutoff temp for H - only material in this temperature range will be "detected"
HI_min = 0      # K
HI_max = 10000  # K
# Assumed distance to the MS - used to calculate angular resolution
MS_dist = 65 # kpc
# Parameters for velocity dispersion calculation
num_bins = 200    # affects emission line resolution. more bins >> higher velocity/frequency resolution
vel_min = -100    # km/s
vel_max = 100     # km/s
bins = np.linspace(vel_min,vel_max,num=num_bins,retstep=True)
################################################################################
def gaussian_curve(x,mu,sigma,A):
    return (A/(sigma*sqrt(2*pi)))*e**(-0.5*((x-mu)/sigma)**2)

def fixed_res_array(file_name,level,order):
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
    for q1 in xrange(all_data_level_x.ActiveDimensions[t1]):
        vbin = []
        for q2 in xrange(all_data_level_x.ActiveDimensions[t2]):
            weights = np.array([])
            weights = np.resize(weights,(num_bins-1,1))
            for q3 in xrange(all_data_level_x.ActiveDimensions[t3]):
                coords = tuple([q1,q2,q3][i] for i in trans_matrix)
                temp = float(all_data_level_x["temperature"][coords])
                if temp <= HI_max and temp >= HI_min:
                    # velocity dispersion calculation
                    cell_weights = np.array([])
                    cell_weights = np.resize(weights,(num_bins-1,1))
                    vel = float(all_data_level_x[LOS_velocity][coords].in_units("km/s"))
                    mass = float(m_avg*all_data_level_x["H_nuclei_density"][coords]*all_data_level_x["cell_volume"][coords])
                    tb = float(sqrt(2*k*temp/m_avg)/1000)

                    mu = vel
                    sigma = tb
                    A = mass

                    for n in xrange(0,len(bins[0])-1):
                        I = quad(gaussian_curve, bins[0][n], bins[0][n+1],args=(mu,sigma,A))
                        cell_weights[n] += I[0]
                        weights[n] += I[0]
                    if (A-sum(weights))/A > 0.05:
                        if vel < (vel_max+vel_min)*0.5:
                            print " Values out of range. Lower velocity limits."
                            print " Out of range vel: " + str(vel) + " km/s"
                        if vel > (vel_max+vel_min)*0.5:
                            print " Values out of range. Raise velocity limits."
                            print " Out of range vel: " + str(vel) + " km/s"
                        sys.exit()

            if sum(weights) == 0: # all cells on sightline have T outside HI_min < T < HI_max
                sigma = 0
            elif sum(weights) != 0:
                mod = GaussianModel()
                pars = mod.guess(weights.flatten(), x=bins[0][:-1])
                out  = mod.fit(weights.flatten(), pars, x=bins[0][:-1])
                result = np.array(out.params)
                avg_vel = result[1]
                diff_sum = 0
                for a,b in zip(bins[0][:-1],weights.flatten()):
                    diff_sum += b*(a-avg_vel)**2
                sigma = sqrt(diff_sum/sum(weights))

            if FWHM == True:
                vbin.append(2*sqrt(2*ln(2))*sigma)
            elif FWHM == False:
                vbin.append(sigma)
        disp_array.append(vbin)
        prog = (q1+1)*100/float(all_data_level_x.ActiveDimensions[t1])
        if prog % 5 < 100/float(all_data_level_x.ActiveDimensions[t1]):
            print "{0:.1f} %".format(prog)

    da = np.array(disp_array)
    print "fixed resolution array created"
    left = ds.domain_left_edge.in_units("kpc")[t2]
    right = ds.domain_right_edge.in_units("kpc")[t2]
    top = ds.domain_left_edge.in_units("kpc")[t1]
    bottom = ds.domain_right_edge.in_units("kpc")[t1]
    domain = [left,right,top,bottom]
    dx = all_data_level_x["dx"]
    return da,domain,dx

def ang_resolution(dist,dx):
    beam_width = (dx[0,0,0]/(dist*3.08568E18))*3600
    if beam_width >= 3600:
        beam_width = beam_width/3600
        print "At MS distance of %d kpc, beam width is %.1f deg." % (dist,beam_width)
    elif beam_width >= 60:
        beam_width = beam_width/60
        print "At MS distance of %d kpc, beam width is %.1f arcmin." % (dist,beam_width)
    else:
        print "At MS distance of %d kpc, beam width is %.1f arcsec." % (dist,beam_width)
    print "Velocity resolution: {0:.1f} km/s".format(bins[1])

def plotdata(data_array,domain_array,fl_nm,level,order):
    matplotlib.rcParams['font.size'] = 10
    fn = random.randint(0,500)
    fig = plt.figure(fn)
    ax = plt.gca()
    im = ax.imshow(data_array,origin = "lower",aspect = "equal",extent=domain_array,cmap="jet",interpolation='none')
    plt.xlabel(order[1] + " (kpc)")
    plt.ylabel(order[0] + " (kpc)")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label("Velocity Dispersion (km/s)",size=7)

    print "plot created. Saving figure..."
    fig_nm = 'velocity_disp_{0}_lvl_{1}.png'.format(str(fl_nm)[-4:],level)
    plt.savefig(fig_nm,dpi=dpi_level,bbox_inches='tight')
    plt.close()
    print "File saved as: " + fig_nm
    return

def save_data(fl_list,level,order):
    # saves sliceplots to separate directory inside working directory
    directory = os.getcwd() + "/vel_disp_lvl_{0}".format(level)
    if not os.path.exists(directory):
        os.makedirs(directory)
    for fl_nm in fl_list:
        source_file = os.getcwd() + "/velocity_{0}disp_{1}_lvl_{2}.png".format(''.join(order[:1]),fl_nm[-4:],level)
        destination = os.getcwd() + "/vel_{0}disp_lvl_{2}/velocity_{0}disp_{1}_lvl_{2}.png".format(''.join(order[:1]),fl_nm[-4:],level)
        if os.path.exists(source_file):
            print "moving {}...".format(source_file)
            os.rename(source_file, destination)
    print "Data saved to: {1}/vel_{0}disp_lvl_{2}".format(''.join(order[:1]),os.getcwd(),level)

def main((fl_nm,level,xyz_order)):
    da,domain,dx = fixed_res_array(fl_nm,level,xyz_order)
    ang_resolution(MS_dist,dx)
    plotdata(da,domain,fl_nm,level,xyz_order)
    return

if __name__ == "__main__":
    full = raw_input("run all files? (y/N): ")
    fl_nm = os.getcwd() +  "/" + raw_input("enter a filename: ").strip()
    level = int(raw_input("resolution level (0-4): ").strip())
    view = raw_input("XY view? (y/N): ")
    if view == "y":
        xyz_order = ['x','y','z']
    else:
        xyz_order = ['x','z','y']
    if full == "y":
        file_list = glob.glob(fl_nm[0:-4] + "*")
    else:
        file_list = [fl_nm]

    arg_list = zip(file_list,[level]*len(file_list),[xyz_order]*len(file_list))
    start_time = time.time()
    cpu_cores = mp.cpu_count()
    pool = mp.Pool(processes=(cpu_cores-1))
    result = pool.map(main,arg_list)
    pool.close()
    pool.join()
    save_data(file_list,level,xyz_order)
    stop_time = time.time()
    print timer(stop_time-start_time)
