# Running in Parallel

This code may be run serially or in parallel. Running in parallel uses the mpi4py module which must be installed if it is not already.
This can be done by running: 

    pip install mpi4py

To run a script in parallel, run 

    mpirun -np 4 python script.py

where 4 is the number of processors to run on and script.py is the script to run.

More info about Running yt in Parallel can be found here: http://yt-project.org/doc/analyzing/parallel_computation.html#parallel-computation

# Issues

There are a couple of known bugs in 3fieldmaps.py you may run into that are in the process of being fixed now.

**1)**
If your dataset files are not sequential (ex. 0,50,100,150 Myr), or if the time of the last file is greater than the number of dataset file you are running (ex. 20-100 Myr), not all files will be saved to the subdirectory when running all files.

**2)**
If you find yourself running a single file in parallel (for some reason), sometimes the program will finish what it's doing but not exit out. This has to do with the config.txt file being open on threads other than the root processor. You can control-C out of it with your sliceplot intact but mpi4py will get mad at you the next time it tries to run something in parallel because the config.txt file wasn't cleared. Simply run the same command again and it should work.
