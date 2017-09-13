import yt, time
from setup import *

start = raw_input("Start time(Myr): ")
end = raw_input("End time(Myr): ")
if end == "":
    end = start

for x in xrange(int(start),int(end)+1):
    t0 = float(time.time())
    fl_nm = "sing_hdf5_plt_cnt_%04d" % float(x)
    ds = yt.load(fl_nm)

    a = (0,0,5.77122257782E+21)
    b = (0,0,6.32085621864E+21)
    ray = ds.ray((0,0,0),(0,0,3.7028162e22))
    point1 = ds.point(a)
    point2 = ds.point(b)

    vel_diff = (point1["velocity_z"]-point2["velocity_z"]).in_units("km/s")
    avg_vel = ((point1["velocity_z"]+point2["velocity_z"])/2.0).in_units("km/s")
    print "vel_diff " + str(vel_diff)
    print "avg_vel " + str(avg_vel)

    print "cell1 mass: " + str(point1["cell_mass"])
    print "cell2 mass: " + str(point2["cell_mass"])
    ke1 = (point1["kinetic_energy"]*point1["cell_volume"]).in_units("N*m")
    ke2 = (point2["kinetic_energy"]*point2["cell_volume"]).in_units("N*m")
    print "kinetic_energy1: " + str(ke1)
    print "kinetic_energy2: " + str(ke2)
    print "cell1 specific_KE: " + str(point1["specific_KE"])
    print "cell2 specific_KE: " + str(point2["specific_KE"])

    rad_fl = "rad_data_%04d.txt" % float(x)
    vel_fl = "vel_data_%04d.txt" % float(x)
    mass_fl = "mass_data_%04d.txt" % float(x)

    with open (rad_fl, "w") as f:
        for x in ray["radius"]:
            x = str(x)[:-3]
            f.write(x)
            f.write(" ")

    with open (vel_fl, "w") as f:
        for x in ray["velocity_Z"]:
            x = str(x)[:-5]
            f.write(x)
            f.write(" ")

    with open (mass_fl, "w") as f:
        for x in ray["cell_mass"]:
            x = str(x)[:-3]
            f.write(x)
            f.write(" ")

    narf = ((ds.domain_right_edge-ds.domain_left_edge)/2.0)+ds.domain_left_edge
    slc = yt.SlicePlot(ds,"y","velocity_Z",origin = "native",center = [0,0,narf[2]])
    slc.annotate_marker(a,coord_system='data',plot_args={'color':'white','s':10})
    slc.annotate_marker(b,coord_system='data',plot_args={'color':'white','s':10})
    slc.annotate_timestamp(time_format="t = {time:.0f} {units}",draw_inset_box=True,corner='upper_left')
    slc.save()

    prof = yt.create_profile(ray, 'radius', 'velocity_Z',
                                units = {'radius': 'kpc'},
                                logs = {'radius': False})

    plot = yt.ProfilePlot.from_profiles(prof)
    plot.set_log("velocity_Z",False)
    plot.save()
    t1 = float(time.time())
    print "%4d DONE" % float(x)
    timer(t1 - t0)
