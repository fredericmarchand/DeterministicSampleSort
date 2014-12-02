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

    int size;
    int sortedDataSize = 0;

    int *data;
    int *localPSample;
    int *allPSamples;
    int *counts;
    int *recCounts;
    int *sendCounts;
    int *recDisplacements;
    int *sendDisplacements;
    int *bucketLocation;
    int *bucketSize;
    int *sortedData;

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
    
    size = n/p;

    data = new int[size + 1];
    sortedData = new int[(2 * size) +1];
    localPSample = new int[p + 1];
    allPSamples  = new int [p * p + 1];
    counts            = new int[p + 1];
    recCounts         = new int[p + 1];
    recDisplacements  = new int[p + 1];
    sendCounts        = new int[p + 1];
    sendDisplacements = new int[p + 1];
    bucketLocation    = new int[p + 1];
    bucketSize        = new int[p + 1];

    while (getline(input, buffer))
    {
        data[i++] = atoi(buffer.c_str());
    }
    
    input.close();

    //Sort locally
    heapsort(data, size);

    //Create p-sample
    for (int i = 0; i < p; ++i) 
    {
        localPSample[i] = data[i * size / p];
        recDisplacements[i] = i * p;
        recCounts[i] = p;
    }

    //Send all p-samples to proc 1
    MPI_Gatherv(localPSample, p, MPI_INT, allPSamples, recCounts, recDisplacements, MPI_INT, 0, MPI_COMM_WORLD);

    //Sort all received samples and compute global p-sample
    if (id == 0)
    {
        heapsort(allPSamples, p * p);
        for (int i = 0; i < p; ++i)
        {
            localPSample[i] = allPSamples[i * p];
        }
    }
    
    //Broadcast global p-sample
    MPI_Bcast(localPSample, p, MPI_INT, 0, MPI_COMM_WORLD);
    
    //Bucket locally according to global p-sample
    int j = 0;
    bucketLocation[0] = 0;
    for (int i = 1; i < p; ++i)
    {
        while (data[j] < localPSample[i])
        {
            ++j;
        }
        bucketLocation[i] = j; 
    }
    
    //Send bucket i to proc i
    for (int i = 0; i < p-1; ++i)
    {
        bucketSize[i] = bucketLocation[i+1] - bucketLocation[i];
    }
    bucketSize[p-1] = size - bucketLocation[p-1];
    
    for (int i = 0; i < p; ++i)
    {
        sendCounts[i] = 1;
        sendDisplacements[i] = i;
        recCounts[i] = 1;
        recDisplacements[i] = i;
    }

    MPI_Alltoallv(bucketSize, sendCounts, sendDisplacements, MPI_INT, counts, recCounts, recDisplacements, MPI_INT, MPI_COMM_WORLD);

    sortedDataSize = 1;
    for (int i = 0; i < p; ++i)
    {
        sendCounts[i] = bucketSize[i];
        sendDisplacements[i] = bucketLocation[i];
        recCounts[i] = counts[i];
        recDisplacements[i] = sortedDataSize-1; 
        sortedDataSize += counts[i];
    }
    sortedDataSize--;

    MPI_Alltoallv(data, sendCounts, sendDisplacements, MPI_INT, sortedData, recCounts, recDisplacements, MPI_INT, MPI_COMM_WORLD);

    //Resort locally
    heapsort(sortedData, sortedDataSize);

    //balance? 

    //Write to file
    ofstream outputFile (outputFilePath); 

    if (!outputFile.is_open())
    {
        cerr << "Could not open file " << outputFilePath << endl;
        return 1;   /// Error
    }

    for (int i = 0; i < sortedDataSize; ++i)
    {
        outputFile << sortedData[i] << endl;
    }

    outputFile.close();


    // Terminate MPI.
    MPI::Finalize();

    delete [] data;
    delete [] localPSample;
    delete [] allPSamples;
    delete [] counts;
    delete [] recCounts;
    delete [] sendCounts;
    delete [] recDisplacements;
    delete [] sendDisplacements;
    delete [] bucketLocation;
    delete [] bucketSize;
    delete [] sortedData;

    return 0;
}
