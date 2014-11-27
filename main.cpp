#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <fstream>
                   
using namespace std;

#include "mpi.h"

#define DEBUG 2 

void swap(int *array1, int *array2)
{
    int temp = array1[0];
    array1[0] = array2[0];
    array2[0] = temp;
}

void siftDown(int array[], int start, int end)
{
    int root = start;
            
    while (root*2+1 < end)
    {
        int child = 2*root + 1;
        if ((child + 1 < end) && (array[child] < array[child+1]))
        {
            child += 1;
        }
        if (array[root] < array[child])
        {
            swap(&array[child], &array[root]);
            root = child;
        }
        else
            return;
    }
}

//Worst Case O(n*log(n))
void heapsort(int array[], int count)
{
    int start;
    int end;
                
    for (start = (count-2)/2; start >= 0; start--)
    {
        siftDown(array, start, count);
    }
    for (end = count-1; end > 0; end--)
    {
        swap(&array[end], &array[0]);
        siftDown(array, 0, end);
    }
}

int main(int argc, char *argv[])
{
    char inputFilePath[64];
    char outputFilePath[64];
    int n;
    int p;
    ifstream input;
    string buffer;
    int i;

    int *inputArray;
    int *recBuffer;

    int processors;
    int id;

    MPI::Init(argc, argv); //  Initialize MPI.
    processors = MPI::COMM_WORLD.Get_size(); //  Get the number of processes.
    id = MPI::COMM_WORLD.Get_rank(); //  Get the individual process ID.

    sprintf(inputFilePath, "input-%d.txt", id);
    sprintf(outputFilePath, "output-%d.txt", id);
    
    input.open(inputFilePath);

    if (!input.is_open())
    {
        cerr << "Could not open file " << inputFilePath << endl;
        return 1;   /// Error
    }
    
    getline(input, buffer);
    n = atoi(buffer.c_str());
    
    getline(input, buffer);
    p = atoi(buffer.c_str());
    
    inputArray = new int[n/p];
    recBuffer = new int[n/p];
    while (getline(input, buffer))
    {
        inputArray[i++] = atoi(buffer.c_str());
    }
    
    input.close();

#if DEBUG == 1
    for (int i = 0; i < n/p; ++i)
        cout << inputArray[i] << ", ";
    cout << endl;
#endif

    heapsort(inputArray, n/p);

#if DEBUG == 1
    for (int i = 0; i < n/p; ++i)
        cout << inputArray[i] << ", ";
    cout << endl;
#endif

//    MPI_Alltoall(inputArray, n/p, MPI_INT, recBuffer, n/p, MPI_INT, MPI_COMM_WORLD);

    if (id == 0){
#if DEBUG == 2
    for (int i = 0; i < n/p; ++i)
        cout << recBuffer[i] << ", ";
    cout << endl;
#endif

    heapsort(recBuffer, n/p);

#if DEBUG == 2 
    for (int i = 0; i < n/p; ++i)
        cout << recBuffer[i] << ", ";
    cout << endl;
#endif
    }




    // Terminate MPI.
    MPI::Finalize();

    delete [] inputArray;
    delete [] recBuffer;

    return 0;
}
