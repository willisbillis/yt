import yt, os, sys
from PIL import Image
from yt import YTQuantity
from setup import *

################################################################################

def multi():
    start = raw_input("Start time(Myr): ")
    end = raw_input("End time(Myr): ")
    if end == "":
        end = start

    # make new directory for images
    directory = os.getcwd() + "/"
    save_dir = directory + "Multiplot/"
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # check if sliceplots exist for fields specified
    fields = [['density'],['velocity_z']]
    for i in fields:
        field_dir = directory + i[0] + "_y/"
        fields[fields.index(i)].append(field_dir)
        if not os.path.exists(field_dir):
            print field_dir + " does not exist"
            sys.exit()

    pfx = "dual_4r_hdf5_plt_cnt_"
    sfx = "_Slice_y_"
    for i in xrange(int(start),int(end)+1):
        num = "%04d" % i
        im = []

        for field in fields:
            slc_path = field[1] + pfx + num + sfx + field[0] + ".png"
            im.append(Image.open(slc_path))

        vert = 265
        horz = 1060

        new_im = Image.new("RGB", (1060,530))

        for i in range(2):
            new_im.paste(im[i], (0,i*vert))

        new_im.save(save_dir + "multi_" + num + ".png")

multi()
