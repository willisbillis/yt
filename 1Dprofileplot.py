import yt, os, sys, time, traceback, threading

yt.enable_parallelism()

# This is an interactive UI for loading in data sets to yt and creating 1D PhasePlots.
# It takes raw_inputs from the user to determine the file(s) and fields to
# run the PhasePlots on. These inputs are written to a configuration file so it
# can be accessed by multiple processors when run in parallel. Every file is
# checked to make sure that the fields are valid inputs by the debug()
# function. Tags such as 'CONFIG OPEN', 'debug', and 'COMPLETE (2.2)' are
# printed to update the user on what is happening in the program.

    #insert into script to print full callstack (for debugging)
def trace():
    for line in traceback.format_stack():
        print line.strip()

    #assigns variable to config file and returns open file
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

    #wipes config file after use in the event the next run would try to use old config parameters
def clear_config():
    config = open("config.txt", "w")
    config.close()

    #write input parameters to a configuration txt file so all processors have access
def input():
    print "input"
    config = open("config.txt", "w")
    print "OPEN (1)"
    full = raw_input("run all files? (y/n): ")
    fl_nm = raw_input("enter filename: ").strip()
    field_x = raw_input("first field: ")
    field_y = raw_input("second field: ")
    if full != "y":
        full = "n"
    config.write(full + '\n')
    config.write(fl_nm + '\n')
    config.write(field_x + '\n')
    config.write(field_y)
    config.close()
    print "CLOSED (1)"

    #checks for given fields in stored yt field lists
def debug(loadfile,field_x,field_y):
    print "debug"
    fields = [x[1] for x in loadfile.field_list] + [x[1] for x in loadfile.derived_field_list]
    if field_x in fields:
        return True
    else:
        print "%s not available (field_x)" % field_x
        return False
    if field2 in fields:
        return True
    else:
        print "%s not available (field_y)" % field_y
        return False

def time_fix(t):
    if t <= 60:
        print "Elapsed time: %.1f sec" % t
    elif t <= 7200:
        t = float(t)/60
        print "Elapsed time: %.1f min" % t
    return

def save_data(field_x,field_y):
    print "save_data"
    directory = "/data/mewilliams/Run1/1DProfile_{0}_{1}".format(field_x,field_y)
    if not os.path.exists(directory):
        os.makedirs(directory)
    source_file = "/data/mewilliams/Run1/Multi-data_1d-Profile_{0}_{1}.png".format(field_x,field_y)
    destination = "/data/mewilliams/Run1/{0}_{1}/Multi-data_1d-Profile_{0}_{1}.png".format(field_x,field_y)
    if os.path.exists(source_file):
        os.rename(source_file, destination)
    print "Data saved to: /data/mewilliams/Run1/1DProfile_{0}_{1}".format(field_x,field_y)

    #rescales x-axis of final plot for consistency
def rescale(plot,field_x):
    if field_x == "density":
        plot.set_xlim(1e-28,1e-22)
    if field_x == "temperature":
        plot.set_xlim(1e2,1e7)
    return plot

    #run all files
def run_set(config_fl):
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
        profiles = []
        labels = []
        for ds in ts.piter():
            ad = ds.all_data()
            if yt.is_root():
                # Create a 1d profile of field_x vs. field_y
                profiles.append(yt.create_profile(ad,field_x,fields=[field_y]))
                # Add labels
                labels.append("t = %.f Myr" % ds.current_time.in_units("Myr"))
                # Create the profile plot from the list of profiles
                plot = yt.ProfilePlot.from_profiles(profiles, labels=labels)
                rescale(plot,field_x)
                # Save the image
                plot.save()
    else:
        print "NOT RUN"
        sys.exit()
    if yt.is_root():
        print "FILES RUN: " + str(len(ts))
        print "FIELD_X: " + field_x
        print "FIELD_Y: " + field_y
        save_data(field_x,field_y,ts)
    t1 = float(time.time())
    if yt.is_root():
        time_fix(t1 - t0)

#run program

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
            trace()
            rescale(plot,field_x)
            plot.save()
        print "FILE NAME: " + fl_nm
        print "FIELD_X: " + field_x
        print "FIELD_Y: " + field_y
        t1 = float(time.time())
        print "COMPLETE (2.2)"
        time_fix(t1 - t0)

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
