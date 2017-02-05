import yt, os, sys, time, traceback, threading

yt.enable_parallelism()

# This is an interactive UI for loading in data sets to yt and performing SlicePlots.
# It takes raw_inputs from the user to determine the file(s), field, and axis to
# run the SlicePlots on. These inputs are written to a configuration file so it
# can be accessed by multiple processors when run in parallel. Every file is
# checked to make sure that the field and axes are valid inputs by the debug()
# function. Tags such as 'CONFIG OPEN', 'debug', and 'COMPLETE (2.2)' are
# printed to update the user on what is happening in the program.

    # insert into script to print full callstack (for debugging)
def trace():
    for line in traceback.format_stack():
        print line.strip()

    # assigns variable to config file and returns open file
def open_config():
    config = open("config.txt", "r")
    config.close()
    if config.closed:
        print "CONFIG OPEN // ON THREAD: " + str(threading.current_thread())
        config = open("config.txt", "r")
        return config
    else:
        print "ERROR: CANNOT READ CONFIG"
        sys.exit()

    # wipes config file after use in the event the next run would try to use old config parameters
def clear_config():
    config = open("config.txt", "w")
    config.close()

    # write input parameters to a configuration txt file so all processors have access
def input():
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

    # checks for given field in stored yt field lists and given axis in axis list
def debug(loadfile,field,axis):
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

    # outputs time in more readable units
def time_fix(t):
    if t <= 60:
        print "Elapsed time: %.1f sec" % t
    elif t <= 7200:
        t = float(t)/60
        print "Elapsed time: %.1f min" % t
    return

    # saves sliceplots to separate directory inside working directory
def save_data(ds,field,axis,tseries):
    directory = os.getcwd() + "/{1}_{0}".format(axis,field)
    prefix = ds[0:-4]
    if not os.path.exists(directory):
        os.makedirs(directory)
    for x in xrange(0,len(tseries)):
        print "checking file %04d" % x
        source_file = os.getcwd() + "/{0}%04d_Slice_{1}_{2}.png".format(prefix,axis,field) % x
        destination = os.getcwd() + "/{2}_{1}/{0}%04d_Slice_{1}_{2}.png".format(prefix,axis,field) % x
        if os.path.exists(source_file):
            print "moving %s..." % source_file
            os.rename(source_file, destination)
    print "Data saved to: %s/{1}_{0}".format(axis,field) % os.getcwd()

    # runs a single file
def runfile(loadfile,field,axis):
    print "runfile"
    narf = ((loadfile.domain_right_edge-loadfile.domain_left_edge)/2.0)+loadfile.domain_left_edge
    slc = yt.SlicePlot(loadfile,axis,[field],origin = "native",center = [0,0,narf[2]])
        # resizes axis of metallicity plot to limits defined in Gritton et al., 2014
    if field == "met_O":
        slc.set_zlim("met_O",0.0001,1.0)
    if yt.is_root():
        slc.annotate_timestamp()
        slc.save()

    # run all files
def run_set(config_fl):
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

# run program

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
