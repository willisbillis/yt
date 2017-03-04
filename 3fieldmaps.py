import yt, os, sys, time, traceback, threading
from yt import YTQuantity
from yt.units import dimensions

yt.enable_parallelism()

# This is an interactive UI for loading in data sets to yt and performing SlicePlots.
# It takes raw_inputs from the user to determine the file(s), field, and axis to
# run the SlicePlots on. These inputs are written to a configuration file so it
# can be accessed by multiple processors when run in parallel. Every file is
# checked to make sure that the field and axes are valid inputs by the debug()
# function. Tags such as 'CONFIG OPEN', 'debug', and 'COMPLETE (2.2)' are
# printed to update the user on what is happening in the program.

################################################################################
    # insert into script to print full callstack (for debugging)
def trace():
    for line in traceback.format_stack():
        print line.strip()

################################################################################
    # new fields
def cloud_velz(field,data):
    bulk_vel=YTQuantity(140e5, 'cm/s')
    return data['flash', 'velz'].in_units('cm/s')-bulk_vel
yt.add_field(("gas", "cloud_velz"), units="cm/s", function=cloud_velz)

def dyn_pressure(field,data):
    dyn_press=data['gas', 'kinetic_energy']/data['gas', 'cell_volume']
    return dyn_press
yt.add_field(("gas", "dyn_pressure"), units="auto", dimensions=dimensions.pressure, function=dyn_pressure)

def partcntH(field,data):
    molmass = YTQuantity(1.00794, "g")
    partcnt = data['flash','h   ']*data['gas','cell_mass']/molmass*6.0221409e23
    return partcnt
yt.add_field(('gas','partcntH'), function = partcntH, units = "")

def partcntO(field,data):
    molmass = YTQuantity(15.9994, "g")
    partcnt = (data['flash','o   '] +
                data['flash','o1  '] +
                data['flash','o2  '] +
                data['flash','o3  '] +
                data['flash','o4  '] +
                data['flash','o5  '] +
                data['flash','o6  '] +
                data['flash','o7  '] +
                data['flash','o8  '] )*data['gas','cell_mass']/molmass*6.0221409e23
    return partcnt
yt.add_field(('gas','partcntO'), function = partcntO, units = "")

def pressure(field,data):
    R = YTQuantity(8.3144598e6, 'cm**3*Pa/K')
    pressure = (data['gas','partcntO']+data['gas','partcntH'])/6.0221409e23*R*data['gas','temperature'].in_units("K")/data['gas','cell_volume'].in_units("cm**3")
    return pressure.convert_to_units("Pa")
yt.add_field(("gas","pressure"), units = "Pa", function=pressure, force_override=True)

def ram_pressure_z(field,data):
    ram_press=0.5*data['gas', 'density']*(data['gas', 'velocity_z']**2)
    return ram_press
yt.add_field(("gas", "ram_pressure_z"), units="g*cm**-1*s**-2", function=ram_pressure_z, force_override=True)

def ram_pressure_tot(field,data):
    ram_press=0.5*data['gas','density']*((data['gas', 'velocity_x']**2)+(data['gas', 'velocity_y']**2)+(data['gas', 'velocity_z']**2))
    return ram_press
yt.add_field(("gas", "ram_pressure_tot"), units="g*cm**-1*s**-2", function=ram_pressure_tot, force_override=True)

def mach_speed(field,data):
    sound_speed = (data['gas','pressure']/data['gas','density'].in_units('kg/m**3'))**0.5
    mach_speed = (((data['gas', 'velocity_x']**2)+(data['gas', 'velocity_y']**2)+(data['gas', 'velocity_z']**2))**(0.5)).in_units("m/s")/sound_speed
    return mach_speed
yt.add_field(('gas','mach_speed'), units="", function=mach_speed)

################################################################################
def open_config():
    # assigns variable to config file and returns open file
    config = open("config.txt", "r")
    config.close()
    if config.closed:
        print "CONFIG OPEN // ON THREAD: " + str(threading.current_thread())
        config = open("config.txt", "r")
        return config
    else:
        print "ERROR: CANNOT READ CONFIG"
        sys.exit()
def clear_config():
    # wipes config file after use in the event the next run would try to use old config parameters
    config = open("config.txt", "w")
    config.close()
def input():
    # write input parameters to a configuration txt file so all processors have access
    print "input"
    config = open("config.txt", "w")
    print "OPEN (1)"
    full = raw_input("run all files? (y/n): ")
    fl_nm = raw_input("enter filename: ").strip()
    field = raw_input("field to run: ")
    axis = raw_input("axis(if y hit return): ")
    if axis == "":
        axis = "y"
    if full != "y":
        full = "n"
    config.write(full + '\n')
    config.write(fl_nm + '\n')
    config.write(field + '\n')
    config.write(axis)
    config.close()
    print "CLOSED (1)"

def debug(loadfile,field,axis):
    #checks for given fields in stored yt field lists
    print "debug"
    fields = [x[1] for x in loadfile.field_list] + [x[1] for x in loadfile.derived_field_list]
    if field in fields:
        axes = ["x","y","z"]
        if axis not in axes:
            print "ERROR: invalid axis"
            return False
        return True
    else:
        print "ERROR: field not available"
        return False
def rescale(plot,field):
    if field == "met_O":
        # resizes axis of metallicity plot to limits defined in Gritton et al., 2014
        plot.set_zlim("met_O",0.0001,1.0)
    if field == "cloud_velz":
        plot.set_zlim("cloud_velz",-1e8,1e7)
    #if field == "density":
    #    plot.set_zlim("density",1e-27,1e-22)
    return plot
def runfile(loadfile,field,axis):
    # runs a single file
    print "runfile"
    narf = ((loadfile.domain_right_edge-loadfile.domain_left_edge)/2.0)+loadfile.domain_left_edge
    slc = yt.SlicePlot(loadfile,axis,[field],origin = "native",center = [0,0,narf[2]])
    rescale(slc,field)
    if yt.is_root():
        slc.annotate_timestamp(time_format="t = {time:.0f} {units}")   # You can customize your slices here (ex. add contours,etc.)
        slc.save()

def save_data(ds,field,axis,tseries):
    # saves sliceplots to separate directory inside working directory
    directory = os.getcwd() + "/{1}_{0}".format(axis,field)
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
def time_fix(t):
    # outputs time in more readable units
    if t <= 60:
        print "Elapsed time: %.1f sec" % t
    elif t <= 7200:
        t = float(t)/60
        print "Elapsed time: %.1f min" % t
    return
def run_set(config_fl):
    # run all files
    print "run_set"
    t0 = float(time.time())
    # first line of config file must be stored as null so that fl_nm,field, and
    # axis are read from the correct lines
    null = config_fl.readline().strip()
    fl_nm = config_fl.readline().strip()
    field = config_fl.readline().strip()
    axis = config_fl.readline().strip()

    # run all files with format dataset_prefix_XXXX
    time_series = fl_nm[0:-4] + "????"
    ts = yt.load(time_series)
    ds = yt.load(fl_nm)
    if debug(ds,field,axis) == True:
        for ds in ts.piter():
            runfile(ds,field,axis)
    else:
        print "       check field/axis inputs ({0},{1})".format(field,axis)
        clear_config()
        sys.exit()
    if yt.is_root():
        save_data(fl_nm,field,axis,ts)
        print "FILES RUN: " + str(len(ts))
        print "FIELD: " + field
        print "AXIS: " + axis
    t1 = float(time.time())
    if yt.is_root():
        time_fix(t1 - t0)
################################################################################
# Running the program

# input function is run serially to only write parameters to config file once
if yt.is_root():
    input()

def main_set(config_fl):
    if yt.is_root():
        print "START (2.1)"
    run_set(config_fl)
    if yt.is_root():
        print "COMPLETE (2.1)"

def main_sing(config_fl):
    if yt.is_root():
        print "START (2.2)"
        t0 = float(time.time())
        null = config_fl.readline().strip()
        fl_nm = config_fl.readline().strip()
        field = config_fl.readline().strip()
        axis = config_fl.readline().strip()
        ds = yt.load(fl_nm)
        if debug(ds,field,axis) == True:
            runfile(ds,field,axis)
        else:
            print "       check field/axis inputs ({0},{1})".format(field,axis)
            clear_config()
            sys.exit()
        print "FIELD: " + field
        print "AXIS: " + axis
        t1 = float(time.time())
        print "COMPLETE (2.2)"
        time_fix(t1 - t0)
        print "OUTPUT FILE NAME: " + fl_nm + "_Slice_{1}_{0}.png".format(field,axis)

while True:
    config = open("config.txt","r")
    full = config.readline().strip()
    config.close()
    if yt.is_root():
        print "FULL: " + str(full)
    full_list = ["y","n"]
    if full not in full_list:
        time.sleep(10)
        if yt.is_root():
            print "waiting"
    elif full == "y":
        parameters = open_config()
        main_set(parameters)
        clear_config()
        break
    elif full == "n":
        parameters = open_config()
        main_sing(parameters)
        clear_config()
        break
