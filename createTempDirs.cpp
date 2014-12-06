#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <fstream>
#include <sys/stat.h>

using namespace std;

#define MASK 0777

#include "mpi.h"

int main(int argc, char *argv[])
{
	char path[256];
	char filename[256];
	char permissions[256];

	int procs;
	int id;

	int n, p;

	MPI::Init(argc, argv); //  Initialize MPI

	if (argc != 3)
	{
		MPI::Finalize();		
		return 1;
	}

	procs = MPI::COMM_WORLD.Get_size(); //  Get the number of processes.
	id = MPI::COMM_WORLD.Get_rank(); //  Get the individual process ID.

	n = atoi(argv[1]);
	p = atoi(argv[2]);

	sprintf(path, "/tmp/fredericmarchand");
	mkdir(path, MASK); 

	sprintf(filename, "%s/input-%d.txt", path, id);

	//Write to file
	ofstream outputFile (filename);
	
	if (!outputFile.is_open())
	{
		cerr << "Could not open file " << filename << endl;
		return 1;   /// Error
	}
	
	outputFile << n << endl;
	outputFile << p << endl;

	for (int i = 0; i < n/p; ++i)
	{
		outputFile << rand() % 100000 << endl;
	}
	
	outputFile.close();

	sprintf(permissions, "chmod 777 %s", path);
	system(permissions);
	sprintf(permissions, "chmod 777 %s", filename);
	system(permissions);

	// Terminate MPI.
	MPI::Finalize();

	return 0;
}
