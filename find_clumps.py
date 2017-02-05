import numpy as np
import time

import yt
from yt import YTQuantity
from yt.analysis_modules.level_sets.api import *

yt.enable_parallelism()
t0 = time.time()
ds = yt.load("4.2.1_sap_hdf5_plt_cnt_0100")
data_source = ds.all_data()

def thermal_energy(field,data):
    k_boltzmann = YTQuantity(1.38064852e-16, "(cm**2/s**2)/K")
    return 3/2*(data['temperature'])*(k_boltzmann)

yt.add_field(('gas','thermal_energy'), function = thermal_energy, units="cm**2/s**2", force_override=True)

def timer(t):
    if t <= 60:
        return "%.2f sec" % t
    elif t <= 7200:
        t = float(t)/60
        return "%.1f min" % t
    elif t <= 129600:
        t = float(t)/3600
        return "%.2f hours" % t

# the field to be used for contouring
field = ("gas", "density")

# This is the multiplicative interval between contours.
step = 2.0

# Now we set some sane min/max values between which we want to find contours.
# This is how we tell the clump finder what to look for -- it won't look for
# contours connected below or above these threshold values.
c_min = 10**np.floor(np.log10(data_source[field]).min()  )
c_max = 10**np.floor(np.log10(data_source[field]).max()+1)

# Now find get our 'base' clump -- this one just covers the whole domain.
master_clump = Clump(data_source, field)

# Add a "validator" to weed out clumps with less than 20 cells.
# As many validators can be added as you want.
master_clump.add_validator("min_cells", 20)


t1 = time.time()
if yt.is_root():
    print "Part 1: %s" % timer(t1 - t0)

# Begin clump finding.
find_clumps(master_clump, c_min, c_max, step)

t2 = time.time()
if yt.is_root():
    print "Part 2: %s" % timer(t2 - t1)

if yt.is_root():
    # Write out the full clump hierarchy.
    write_clump_index(master_clump, 0, "%s_clump_hierarchy.txt" % ds)

    # Write out only the leaf nodes of the hierarchy.
    write_clumps(master_clump,0, "%s_clumps.txt" % ds)

    # We can traverse the clump hierarchy to get a list of all of the 'leaf' clumps
    leaf_clumps = get_lowest_clumps(master_clump)

    # If you'd like to visualize these clumps, a list of clumps can be supplied to
    # the "clumps" callback on a plot.  First, we create a projection plot:
    narf = ((ds.domain_right_edge-ds.domain_left_edge)/2.0)+ds.domain_left_edge
    prj = yt.ProjectionPlot(ds,"y", field,origin = "native",center = [0,0,narf[2]])

    # Next we annotate the plot with contours on the borders of the clumps
    prj.annotate_clumps(leaf_clumps)

    # Lastly, we write the plot to disk.
    prj.save('clumps')

t3 = time.time()
if yt.is_root():
    print "Part 3: %s" % timer(t3 - t2)
    print "Total time: %s" % timer(t3 - t0)
