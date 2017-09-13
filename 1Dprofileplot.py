import yt, os, sys, time, traceback, threading
from yt import YTQuantity
from yt.units import dimensions
from setup import *

yt.enable_parallelism()

# This is an interactive UI for loading in data sets to yt and creating 1D ProfilePlots.
# It takes raw_inputs from the user to determine the file(s) and fields to
# run the PhasePlots on. These inputs are written to a configuration file so it
# can be accessed by multiple processors when run in parallel. Every file is
# checked to make sure that the fields are valid inputs by the debug()
# function. Tags such as 'CONFIG OPEN', 'debug', and 'COMPLETE (2.2)' are
# printed to update the user on what is happening in the program.

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
    # write input parameters to a configuration txt file so all processors have access in parallel
    print "input"
    config = open("config.txt", "w")
    print "OPEN (1)"
    full = raw_input("run all files? (y/n): ")
    fl_nm = raw_input("enter filename: ").strip()
    field_x = raw_input("x-axis field: ")
    field_y = raw_input("y-axis field: ")
    if full != "y":
        full = "n"
    config.write(full + '\n')
    config.write(fl_nm + '\n')
    config.write(field_x + '\n')
    config.write(field_y)
    config.close()
    print "CLOSED (1)"
def debug(loadfile,field_x,field_y):
    #checks for given fields in stored yt field lists
    print "debug"
    fields = [x[1] for x in loadfile.field_list] + [x[1] for x in loadfile.derived_field_list]
    if field_x in fields:
        return True
    else:
        print "%s not available (field_x)" % field_x
        return False
    if field_y in fields:
        return True
    else:
        print "%s not available (field_y)" % field_y
        return False
def rescale(plot,field_x):
    #rescales x-axis of final plot for consistency
    if field_x == "density":
        plot.set_xlim(1e-28,1e-22)
    if field_x == "temperature":
        plot.set_xlim(1e2,1e7)
    return plot
def runfile(domain_data,field_x,field_y):
    # runs a single file
    print "runfile"
    plot = yt.ProfilePlot(domain_data,field_x,[field_y])
    rescale(plot,field_x)
    if yt.is_root():
        plot.save()
def save_data(ds,field_x,field_y,tseries):
    # saves sliceplots to separate directory inside working directory
    directory = os.getcwd() + "/1D_{0}_{1}".format(field_x,field_y)
    prefix = ds[0:-4]
    if not os.path.exists(directory):
        os.makedirs(directory)
    for x in xrange(0,len(tseries)):
        source_file = os.getcwd() + "/{0}%04d_1d-Profile_{1}_{2}.png".format(prefix,field_x,field_y) % x
        destination = os.getcwd() + "/1D_{1}_{2}/{0}%04d_1d-Profile_{1}_{2}.png".format(prefix,field_x,field_y) % x
        if os.path.exists(source_file):
            print "moving %s..." % source_file
            os.rename(source_file, destination)
    print "Data saved to: %s/1D_{0}_{1}".format(field_x,field_y) % os.getcwd()
def run_set(config_fl):
    #run all files
    print "run_set"
    t0 = float(time.time())
    #first line of config file must be stored as null so that fl_nm,field_x, and
    # field_y are read from the correct lines
    null = config_fl.readline().strip()
    fl_nm = config_fl.readline().strip()
    field_x = config_fl.readline().strip()
    field_y = config_fl.readline().strip()

    data = yt.load(fl_nm)
    #run all files with format 4.2.1_sap_hdf5_plt_cnt_XXXX
    time_series = fl_nm[0:-4] + "????"
    ts = yt.load(time_series)
    if debug(data,field_x,field_y) == True:
        for ds in ts.piter():
            ad = ds.all_data()
            runfile(ad,field_x,field_y)
    else:
        clear_config()
        sys.exit()
    if yt.is_root():
        save_data(fl_nm,field_x,field_y,ts)
        print "FILES RUN: " + str(len(ts))
        print "FIELD_X: " + field_x
        print "FIELD_Y: " + field_y
    t1 = float(time.time())
    if yt.is_root():
        timer(t1 - t0)
################################################################################
#Running the program

#input function is run serially to only write parameters to config file once
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
        field_x = config_fl.readline().strip()
        field_y = config_fl.readline().strip()
        ds = yt.load(fl_nm)
        if debug(ds,field_x,field_y) == True:
            ad = ds.all_data()
            plot = yt.ProfilePlot(ad,field_x,field_y)
            rescale(plot,field_x)
            plot.save()
        print "FILE NAME: " + fl_nm
        print "FIELD_X: " + field_x
        print "FIELD_Y: " + field_y
        t1 = float(time.time())
        print "COMPLETE (2.2)"
        timer(t1 - t0)

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
        if yt.is_root():
            parameters = open_config()
            main_sing(parameters)
            clear_config()
            break
