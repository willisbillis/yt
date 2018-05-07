import yt, time, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from setup import *
from yt.units import kpc


fl = raw_input("Start file: ")
start = fl[-4:]
end = raw_input("End time(Myr): ")
field = raw_input("Field: ")
if end == "":
    end = start

t0 = float(time.time())
for x in xrange(int(start),int(end)+1):
    fl_nm = fl[:-4]+"%04d" % float(x)
    ds = yt.load(fl_nm)

    ray = ds.ray((0,0,ds.domain_left_edge[0]),(0,0,ds.domain_right_edge[2]))

    radius = [float(i) for i in ray["radius"].in_units("kpc")]
    field_data = [float(i) for i in ray[field]]

    plt.plot(radius,field_data)
    plt.xlabel("Radius (kpc)")
    plt.ylabel(field)
    print "Saving plot "+fl_nm+"_1D_profile_"+field+".png"
    plt.savefig(fl_nm+"_1D_profile_"+field+".png")
    plt.cla()

"""
directory = os.getcwd() + "/1D_{0}".format(field)
prefix = ds[0:-4]
if not os.path.exists(directory):
    os.makedirs(directory)
for x in xrange(0,len(tseries)):
    source_file = os.getcwd() + "/{0}%04d_Slice_{1}_{2}.png".format(prefix,axis,field) % x
    destination = os.getcwd() + "/{2}_{1}/{0}%04d_Slice_{1}_{2}.png".format(prefix,axis,field) % x
    if os.path.exists(source_file):
        print "moving %s..." % source_file
        os.rename(source_file, destination)
print "Data saved to: %s/{1}_{0}".format(axis,field) % os.getcwd()
"""
t1 = float(time.time())
print timer(t1 - t0)
