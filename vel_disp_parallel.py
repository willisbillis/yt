from __future__ import division
import yt, os, sys, time, threading
from yt import YTQuantity
from yt.units import dimensions
from setup import *

#yt.enable_parallelism()

# This is an interactive UI for loading in data sets to yt to making velocity
# dispersion plots. It takes raw_inputs from the user to determine the file(s)
# to run. These inputs are written to a configuration file so it can be accessed
# by multiple processors when run in parallel. Tags such as 'CONFIG OPEN' and
# 'COMPLETE (2.2)' are printed to update the user on what is happening in the program.

# set resolution level (dpi)
quick = 150
pro = 750
dpi_level = quick
# Turn on/off FWHM (assumes vel disp along LOS is Gaussian, which is valid)
FWHM = False
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
    full = raw_input("run all files? (y/N): ")
    fl_nm = raw_input("enter filename: ").strip()
    level = raw_input("resolution level: ").strip()
    if full != "y":
        full = "n"
    config.write(full + '\n')
    config.write(fl_nm + '\n')
    config.write(level)
    config.close()
    print "CLOSED (1)"
def runfile(fl_nm,level):
    # runs a single file
    print "runfile"
    da,domain = fixed_res_array(fl_nm,level)
    main(da,domain,fl_nm,level)
    return

def save_data(tseries,level):
    # saves plots to separate directory inside working directory
    directory = os.getcwd() + "/vel_disp_lvl_{0}".format(level)
    if not os.path.exists(directory):
        os.makedirs(directory)
    for x in xrange(0,len(tseries)):
        source_file = os.getcwd() + "/velocity_disp_%04d_lvl_{0}.png".format(level) % x
        destination = os.getcwd() + "/vel_disp_lvl_{0}/velocity_disp_%04d_lvl_{0}.png".format(level) % x
        if os.path.exists(source_file):
            print "moving %s..." % source_file
            os.rename(source_file, destination)
    print "Data saved to: {0}/vel_disp_lvl_{1}".format(os.getcwd(),level)

def save_data_xy(tseries,level):
    # saves plots to separate directory inside working directory
    directory = os.getcwd() + "/vel_xydisp_lvl_{0}".format(level)
    if not os.path.exists(directory):
        os.makedirs(directory)
    for x in xrange(0,len(tseries)):
        source_file = os.getcwd() + "/velocity_xydisp_%04d_lvl_{0}.png".format(level) % x
        destination = os.getcwd() + "/vel_xydisp_lvl_{0}/velocity_xydisp_%04d_lvl_{0}.png".format(level) % x
        if os.path.exists(source_file):
            print "moving %s..." % source_file
            os.rename(source_file, destination)
    print "Data saved to: {0}/vel_xydisp_lvl_{1}".format(os.getcwd(),level)

def run_set(config_fl):
    #run all files
    print "run_set"
    t0 = float(time.time())
    # first line of config file must be stored as null so that fl_nm and level
    # are read from the correct lines
    null = config_fl.readline().strip()
    fl_nm = config_fl.readline().strip()
    level = int(config_fl.readline().strip())

    #run all files with format fl_nm_prefix_XXXX
    time_series = fl_nm[0:-4] + "????"
    ts = yt.load(time_series)
    for ds in ts.piter():
        runfile(ds,level)

    if yt.is_root():
        if xyview == "n":
            save_data(ts,level)
        elif xyview == "y":
            save_data_xy(ts,level)
        print "FILES RUN: " + str(len(ts))
    t1 = float(time.time())
    if yt.is_root():
        timer(t1 - t0)

################################################################################
#Running the program

#input function is run serially to only write parameters to config file once
if yt.is_root():
    input()

xyview = raw_input("view from xy-plane? (y/N): ")
if xyview != "y":
    xyview = "n"
if xyview == "n":
    from vel_disp import *
elif xyview == "y":
    from vel_disp_xy import *    

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
        level = int(config_fl.readline().strip())
        ds = yt.load(fl_nm)
        da,domain = fixed_res_array(ds,level)
        main(da,domain,fl_nm,level)
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
