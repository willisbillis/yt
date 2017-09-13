import yt,time,multiprocessing
import numpy as np
from multiprocessing import Pool, Process, Array

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
# Benchmarks (single processor):
# Level 1 - 7 sec
# Level 2 - 43 sec
# Level 3 - 5.6 min
################################################################################
# set resolution level (dpi)
quick = 150
pro = 750
dpi_level = quick
# Turn on/off FWHM (assumes vel disp is Gaussian - which is valid)
FWHM = False

def fixed_res_array(ds,level):
    all_data_level_x = ds.covering_grid(level=level,left_edge=ds.domain_left_edge,dims=ds.domain_dimensions*2**level)
    disp_array = []
    ds.periodicity = (True,True,True)
    for x in xrange(0,ds.domain_dimensions[0]*2**level):
        vbin = []
        for z in xrange(0,ds.domain_dimensions[2]*2**level):
            v = []
            m = []
            for y in xrange(0,ds.domain_dimensions[1]*2**level):
                vel = all_data_level_x["velocity_magnitude"][x,y,z].in_units("km/s")
                v.append(vel)
                mass = all_data_level_x["cell_mass"][x,y,z]
                m.append(mass)
            sigma = np.sqrt(np.sum((v - np.average(v,weights=m))**2) / np.size(v))
            if FWHM == True:
                vbin.append(2*sqrt(2*ln(2))*sigma)
            else:
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
    MS_dist = 80 # kpc
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
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
    matplotlib.rcParams['font.family'] = "sans-serif"
    matplotlib.rcParams['font.size'] = 10
    im = plt.imshow(data_array,origin = "lower",aspect = "equal",extent=domain_array)
    plt.xlabel("z (kpc)")
    plt.ylabel("x (kpc)")
    cbar = plt.colorbar(pad=0,shrink=0.2535,aspect=7)
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
    start_time = time.time()
    fl_nm = raw_input("enter filename: ").strip()
    level = int(raw_input("resolution level: ").strip())
    ds = yt.load(fl_nm)
    da,domain = fixed_res_array(ds,level)
    main(da,domain,fl_nm,level)
    stop_time = time.time()
    print timer(stop_time-start_time)
