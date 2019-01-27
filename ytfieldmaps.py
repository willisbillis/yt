"""
Create yt SlicePlots of FLASH simulation data files. User specifies file(s),
field(s), and axis either interactively on command line or in code below. If
multiple cores are available, job will be split amongst them to increase efficiency.

Written by Elliott Williams https://github.com/willisbillis/yt-github
"""

import os
import sys
import time
import glob
from multiprocessing import Pool, cpu_count
import yt

################################################################################
## CONSTANTS ##
INTERACTIVE_MODE = True # set constants here that you don't want to specify
                        # from the command line each time
# set usr_constants variables if interactive mode is off
FULL = "y"
FL_NM = ""
FIELD_LIST = ["density"]
AXIS = "Y"
USER_CONSTANTS = [INTERACTIVE_MODE, FULL, FL_NM, FIELD_LIST, AXIS]
#############

## File input/output
QUICK = 150
PRO = 1250
DPI_LEVEL = QUICK        # resolution level of output png file (dots per inch)

## Parallel processing parameters
CHUNKSIZE = 3            # best if a multiple of processes for efficiency
MAXTASKS = CHUNKSIZE * 5 # max tasks per child before printing results and
                         #      starting new child process
################################################################################
## FUNCTIONS ##
def inputs(usr_constants):
    """
    Used to specify variables for each run. If Interactive Mode is True, then
    variables are specified from command line prompts.
    """
    int_mode = usr_constants[0]

    if int_mode is True:
        full = raw_input("run all files? (y/N): ")
        fl_nm = os.getcwd() +  "/" + raw_input("enter a filename: ").strip()
        field_list = raw_input("field(s) to run: ").split(" ")
        axis = raw_input("slice axis (if y hit return): ")
    else:
        full = usr_constants[1]
        fl_nm = usr_constants[2]
        field_list = usr_constants[3]
        axis = usr_constants[4]

    if full != "y":
        full = "n"
        file_list = [fl_nm]
    elif full == "y":
        file_list = glob.glob(fl_nm[0:-4] + "*[0-9]")
    if axis == "":
        axis = "y"

    return full, fl_nm, field_list, axis, file_list
def timer(t_sec):
    """ Output time in more readable units. """
    if t_sec <= 60:
        readout = "Elapsed time: {:.1f} sec".format(t_sec)
    elif t_sec <= 7200:
        t_min = t_sec / 60.0
        readout = "Elapsed time: {:.1f} min".format(t_min)
    return readout
def rescale(plot, field):
    """ resize axis limits of certain fields """
    if field == "met_O":
        # resizes axis of metallicity plot to limits defined in Gritton et al, 2014
        plot.set_zlim("met_O", 0.0001, 1.0)
    if field == "velocity_Z":
        plot.set_log("velocity_Z", False)
    if field == "mach_speed":
        plot.set_log("mach_speed", False)
    if field == "velocity_z":
        plot.set_zlim("velocity_z", -1.5e7, 1.5e7)
        plot.set_log("velocity_z", False)
    if field == "ram_pressure":
        plot.set_zlim("ram_pressure", 3e-17, 1e-12)
    #if field == "density":
    #    plot.set_zlim("density",1e-27,1e-22)
    return plot
def runfile((file_name, field, axis)):
    """ Create sliceplot for a single file """
    print "runfile"
    loadfile = yt.load(file_name)
    narf = ((loadfile.domain_right_edge-loadfile.domain_left_edge)/2.0)+loadfile.domain_left_edge
    slc = yt.SlicePlot(loadfile, axis, [field], origin="native",
                       center=[0, 0, narf[2]], fontsize=14)
    rescale(slc, field)
    slc.annotate_timestamp(time_format="t = {time:.0f} {units}",
                           draw_inset_box=True, corner='upper_left')
    # You can customize your slices here (ex. add contours,etc.)
    #slc.annotate_streamlines('velocity_z','velocity_x')
    slc.save(mpl_kwargs={"dpi":DPI_LEVEL})
    return
def debug(file_name, field_list, axis):
    """ Checks for given fields in stored yt field lists and if axis is valid. """
    print "debug"
    axes = ["x", "y", "z"]
    if axis not in axes:
        print "ERROR: invalid axis ({0})".format(axis)
        return False
    loadfile = yt.load(file_name)
    yt_fields = [x[1] for x in loadfile.field_list] + [x[1] for x in loadfile.derived_field_list]
    bad_fields = []
    for field in field_list:
        if field in yt_fields:
            continue
        else:
            bad_fields.append(field)
    if bad_fields:
        print "ERROR: following field(s) not available {0}".format(tuple(bad_fields))
        return False
    return True
def save_data(field_list, axis, fl_list):
    """ Saves sliceplots to separate directory inside working directory """
    for field in field_list:
        directory = os.getcwd() + "/{1}_{0}".format(axis, field)
        if not os.path.exists(directory):
            os.makedirs(directory)
        for fl_nm in fl_list:
            source_file = os.getcwd() + "/{0}_Slice_{1}_{2}.png".format(fl_nm, axis, field)
            destination = os.getcwd() + "/{2}_{1}/{0}_Slice_{1}_{2}.png".format(fl_nm, axis, field)
            if os.path.exists(source_file):
                print "moving {}...".format(source_file)
                os.rename(source_file, destination)
        print "Data saved to: {2}/{1}_{0}".format(axis, field, os.getcwd())
################################################################################
## PROCESSING ##
if __name__ == '__main__':
    FULL, FL_NM, FIELD_LIST, AXIS, FILE_LIST = inputs(USER_CONSTANTS)
    START_TIME = time.time()
    CPU_CORES = cpu_count()
    POOL = Pool(processes=(CPU_CORES-1), maxtasksperchild=MAXTASKS)
    if debug(FL_NM, FIELD_LIST, AXIS) is False:
        sys.exit()

    ARG_LIST = []
    for FIELD in FIELD_LIST:
        ARG_LIST.append(zip(FILE_LIST, [FIELD]*len(FILE_LIST), [AXIS]*len(FILE_LIST)))
    ARG_LIST = [item for sublist in ARG_LIST for item in sublist]
    POOL.map(runfile, ARG_LIST, chunksize=CHUNKSIZE)
    if FULL == "y":
        save_data(FIELD_LIST, AXIS, FILE_LIST)

    POOL.close()
    POOL.join()
    STOP_TIME = time.time()
    print "FIELD(S):", FIELD_LIST
    print "AXIS:", AXIS
    print timer(STOP_TIME-START_TIME)
