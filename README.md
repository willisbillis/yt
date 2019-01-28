# Overview

In this repository are tools for post-processing and analysis of FLASH data files. Many of these utilize already existing tools from the post-processing software yt (www.yt-project.org), and automate tasks so that you can create density plots for all of your simulation files at once and have them in their own directory, for example. More computationally intensive scripts have been designed to run in parallel for efficiency, so that they can take advantage of the other cores in the machine and drastically increase the speed of completion.

**3fieldmaps.py**

This script creates yt slice plots for any specified field. It is run from inside the directory of the FLASH files you wish to process with the command

    python 3fieldmaps.py
    
and will give prompts so that you can customize its use. For multiple fields, separate the fields by spaces on the command line. To create higher resolution plots (for a poster or talk), set dpi_level to pro at the beginning of 3fieldmaps.py.

**vel_disp.py**

This script creates velocity dispersion maps of FLASH simulation files using yt. Run it from inside the directory of the files you want to process with the command

    python vel_disp.py
    
There will be prompts you answer at the command line to The time it takes to create the maps increases rapidly, so start with a low resolution (0,1) and move up to higher (2,3) as you need. There are also settings you can change inside the file (like output file quality, velocity resolution).

**fitsvel_disp.py**

# Requirements

The scripts in this repository have the following dependencies

    numpy==1.15.4
    matplotlib==2.2.3
    scipy==1.2.0
    lmfit==0.9.12
    multiprocessing==2.6.2.1
    yt==3.5.0
    astropy==2.0.9

# Issues

There are a couple of known bugs in 3fieldmaps.py and 1Dprofileplot.py you may run into that are in the process of being fixed now.

**1)**
If your dataset files are not sequential (ex. 0,50,100,150 Myr), or if the time of the last file is greater than the number of dataset file you are running (ex. 20-100 Myr), not all files will be saved to the subdirectory when running all files.
