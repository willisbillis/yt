import yt,sys,matplotlib
from scipy.integrate import quad
from math import pi,sqrt,e
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lmfit.models import GaussianModel

# set resolution level of png (dpi)
quick = 150
pro = 1250
dpi_level = pro
# Ionization cutoff temp for H
HI_min = 0
HI_max = 10000
# Parameters for velocity dispersion calculation
num_bins = 300     # affects line resolution. more bins >> higher resolution
vel_min = -100        # km/s
vel_max = 200      # km/s
bins = np.linspace(vel_min,vel_max,num=num_bins,retstep=True)
################################################################################
ds = yt.load("dual_4r_hdf5_plt_cnt_0160")
ray = ds.ray((0,1e20,ds.domain_left_edge[0]),(0,1e20,ds.domain_right_edge[2]))

k = 1.380649e-23 # J/K
m_avg = 1.6737236e-27 # kg (mass of hydrogen atom)
def gaussian_curve(x,mu,sigma,A):
    return (A/(sigma*sqrt(2*pi)))*e**(-0.5*((x-mu)/sigma)**2)
################################################################################

def line_emission(ray):
    weights = np.array([])
    weights = np.resize(weights,(num_bins-1,1))
    cl_cnt = 0
    print len(ray["velocity_z"])
    for i in xrange(0,len(ray["velocity_z"])):
        if i < 300 and i % 50 == 0:
            print i
        if i > 300 and i % 25 == 0:
            print i
        # start loop for every cell here
        temp = float(ray["temperature"][i])
        if temp <= HI_max and temp >= HI_min:
            cl_cnt += 1
            cell_weights = np.array([])
            cell_weights = np.resize(cell_weights,(num_bins-1,1))
            vel = float(ray["velocity_z"][i].in_units("km/s"))
            mass = float(m_avg*ray["H_nuclei_density"][i]*ray["cell_volume"][i])

            mu = vel  # vel
            sigma = sqrt(2*k*temp/m_avg)/1000 # tb (in km/s)
            A = mass   # mass

            for n in xrange(0,len(bins[0])-1):
                I = quad(gaussian_curve, bins[0][n], bins[0][n+1],args=(mu,sigma,A))
                cell_weights[n] += I[0]
                weights[n] += I[0]
            if (A-sum(weights))/A > 0.003:
                print " Values out of range. Expand velocity limits."
                sys.exit()
            #fig = plt.hist(bins[0][:-1],bins=bins[0],weights=cell_weights.flatten())
        else:
            pass

    if sum(weights) == 0: # all cells on sightline have T outside HI_min < T < HI_max
        sigma = 0
    elif sum(weights) != 0:
        mod = GaussianModel()
        pars = mod.guess(weights.flatten(), x=bins[0][:-1])
        print type(weights.flatten())
        out  = mod.fit(weights.flatten(), pars, x=bins[0][:-1])
        #print(out.fit_report(min_correl=0.25))
        result = np.array(out.params)
        mod_sigma = result[0]
        print mod_sigma

        avg_vel = result[1]
        diff_sum = 0
        for a,b in zip(bins[0][:-1],weights.flatten()):
            diff_sum += b*(a-avg_vel)**2
        sigma = sqrt(diff_sum/sum(weights))
        print sigma

    percent_diff = 200*abs(sigma-mod_sigma)/(sigma+mod_sigma)
    print str(percent_diff) + " %"

    print "Velocity resolution: {0:.1f} km/s".format(bins[1])
    print "cells measured: " + str(cl_cnt) + "/" + str(len(ray["velocity_z"]))

    plt.plot(bins[0][:-1],weights.flatten(), label="Emission line")
    plt.plot(bins[0][:-1],out.best_fit, label="Gaussian fit", alpha=0.6)
    plt.xlabel("Velocity (km/s)")
    plt.ylabel("Intensity")
    plt.legend()
    plt.savefig("hist.png", dpi=dpi_level)

line_emission(ray)
