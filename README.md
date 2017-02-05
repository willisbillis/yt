# yt

Running in Parallel

This code may be run serially or in parallel. Running in parallel uses the mpi4py module which must be installed if it is not already.
This can be done by running: 

pip install mpi4py

To run a script in parallel, run "mpirun -np 4 python script.py", where 4 is the number of processors to run on and script.py is the script to run.

More info about Running yt in Parallel can be found here: http://yt-project.org/doc/analyzing/parallel_computation.html#parallel-computation

# Issues

There are a couple of known bugs in 3fieldmaps.py you may run into that are in the process of being fixed now.





