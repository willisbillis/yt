import yt,time,multiprocessing,math
import numpy as np
from multiprocessing import Pool, Process, Array
from yt.units import kilometer,second

def timer(t):
    # outputs time in more readable units
    if t <= 60:
        return "Elapsed time: %.1f sec" % t
    elif t <= 7200:
        t = float(t)/60
        return "Elapsed time: %.1f min" % t

# NOTE you can definitely parallelize this for higher resolution if you split up
# the task of integrating and calculating the velocity dispersion for each (x,z)
# coordinate across multiple instances/threads/processors
#
# Benchmarks (single processor of Intel Q9450 @ 2.66GHz):
# Level 1 - 7 sec
# Level 2 - 43 sec
# Level 3 - 5.6 min
################################################################################
# set resolution level (dpi)
quick = 150
pro = 750
dpi_level = quick
# Turn on/off FWHM (assumes vel disp along LOS is Gaussian, which is valid)
FWHM = False
# Ionization cutoff temp for H
HI_min = 6000
HI_max = 20000

def fixed_res_array(ds,level):
    all_data_level_x = ds.covering_grid(level=level,left_edge=ds.domain_left_edge,dims=ds.domain_dimensions*2**level)
    disp_array = []
    ds.periodicity = (True,True,True)
    k = 1.380649e-23 # J/K
    m_avg = 1.6737236e-27 # kg (mass of hydrogen atom)
    for x in xrange(0,ds.domain_dimensions[0]*2**level):
        vbin = []
        for z in xrange(0,ds.domain_dimensions[2]*2**level):
            v = []
            m = []
            for y in xrange(0,ds.domain_dimensions[1]*2**level):
                temp = all_data_level_x["temperature"][x,y,z]
                if temp <= HI_max and temp >= HI_min:
                    #print [x,y,z]
                    vel = all_data_level_x["velocity_y"][x,y,z].in_units("km/s")
                    mass = all_data_level_x["cell_mass"][x,y,z]
                    tb = math.sqrt(2*k*temp/m_avg)/1000*kilometer/second
                    v.append(vel+tb)
                    m.append(mass)


            if sum(m) == 0: # all cells on sightline have T outside HI_min < T < HI_max
                sigma = 0
            elif sum(m) != 0:
                sigma = np.sqrt(np.sum((v - np.average(v,weights=m))**2) / np.size(v))
            if FWHM == True:
                vbin.append(2*sqrt(2*ln(2))*sigma)
            elif FWHM == False:
                vbin.append(sigma)
            disp_array.append(vbin)
        print "{0:.1f} %".format((x+1)*100/float(16*2**level))

    da = np.array(disp_array)
    print "fixed resolution array created"
    left = ds.domain_left_edge.in_units("kpc")[2]
    right = ds.domain_right_edge.in_units("kpc")[2]
    top = ds.domain_left_edge.in_units("kpc")[0]
    bottom = ds.domain_right_edge.in_units("kpc")[0]
    domain = [left,right,top,bottom]
    MS_dist = 65 # kpc
    beam_width = (all_data_level_x["dx"][0,0,0]/(MS_dist*3.08568E18))*3600
    if beam_width >= 3600:
        beam_width = beam_width/3600
        print "At MS distance of %d kpc, beam width is %.1f deg." % (MS_dist,beam_width)
    elif beam_width >= 60:
        beam_width = beam_width/60
        print "At MS distance of %d kpc, beam width is %.1f arcmin." % (MS_dist,beam_width)
    else:
        print "At MS distance of %d kpc, beam width is %.1f arcsec." % (MS_dist,beam_width)
    return da,domain

def main(data_array,domain_array,fl_nm,level):
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    mpl.rcParams['font.size'] = 10
    im = plt.imshow(data_array,origin = "lower",aspect = "equal",extent=domain_array,cmap="jet")
    plt.xlabel("z (kpc)")
    plt.ylabel("x (kpc)")
    cbar = plt.colorbar(pad=0,shrink=0.2275,aspect=6)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label("Velocity Dispersion (km/s)",size=7)

    print "plot created. Saving figure..."
    fig_nm = 'velocity_disp_{0}_lvl_{1}.png'.format(str(fl_nm)[-4:],level)
    plt.savefig(fig_nm,dpi=dpi_level,bbox_inches='tight')
    plt.close()
    print "File saved as: " + fig_nm
    return

"""
if __name__ == "__main__":
    pool = multiprocessing.Pool(4)
    pool.map(main,da)
"""

if __name__ == "__main__":
    fl_nm = raw_input("enter filename: ").strip()
    level = int(raw_input("resolution level: ").strip())
    start_time = time.time()
    ds = yt.load(fl_nm)
    da,domain = fixed_res_array(ds,level)
    main(da,domain,fl_nm,level)
    stop_time = time.time()
    print timer(stop_time-start_time)
