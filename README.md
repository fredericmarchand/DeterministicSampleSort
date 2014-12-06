Frederic Marchand
100817579
Comp4009 A4
Deterministic Sample Sort (MPI)

To compile, run 

> make all

To generate the number files run:

> mpirun -np < p > --hostfile hostfile createDirs < n > < p >  

That will create folders in /tmp/fredericmarchand for each processor and create files named input-i.txt (as per the assignment description)

To run the program run:

> mpirun -np < p > --hostfile hostfile main


//Note: all < p >'s must be the same
