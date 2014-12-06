#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <fstream>
#include <sys/stat.h>

using namespace std;

#include "mpi.h"

int main(int argc, char *argv[])
{
    char path[256];
    ifstream input;
    string buffer;

    int p;
    int id;

    MPI::Init(argc, argv); //  Initialize MPI.
    p = MPI::COMM_WORLD.Get_size(); //  Get the number of processes.
    id = MPI::COMM_WORLD.Get_rank(); //  Get the individual process ID.

    sprintf(path, "/tmp/fredericmarchand");
    mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); 
    
    // Terminate MPI.
    MPI::Finalize();

    return 0;
}
