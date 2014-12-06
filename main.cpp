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
        {
            return;
        }
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

void setupCountsAndDisplacements(int *sCounts, int *rCounts, int *sDispl, int *rDispl, int sc, int rc, int sd, int rd, int p)
{
    for (int i = 0; i < p; ++i)
    {
        sCounts[i] = sc;

        if (sd < 0)
            sDispl[i] = i;
        else
            sDispl[i] = sd;

        rCounts[i] = rc;

        if (rd < 0)
            rDispl[i] = i;
        else
            rDispl[i] = rd;
    }
}

void deterministicSampleSort(int *data, int size, int *finalData, int p, int id, double *compTime, double *commTime)
{
    int sortedDataSize = 0;

    int *localPSample;
    int *allPSamples;
    int *allBuckets;
    int *recCounts;
    int *sendCounts;
    int *recDisplacements;
    int *sendDisplacements;
    int *bucketLocation;
    int *bucketSize;
    int *sortedData;
    int *balancingData;
    int *balancedData;

    double compTime2;
    double commTime2;

    compTime2 = MPI::Wtime();

    sortedData = new int[(2 * size) + 1];
    localPSample = new int[p + 1];
    allPSamples  = new int[p * p + 1];
    allBuckets        = new int[p + 1];
    recCounts         = new int[p + 1];
    recDisplacements  = new int[p + 1];
    sendCounts        = new int[p + 1];
    sendDisplacements = new int[p + 1];
    bucketLocation    = new int[p + 1];
    bucketSize        = new int[p + 1];
    balancingData     = new int[p + 1];
    balancedData      = new int[p + 1];
    
    //Sort locally
    heapsort(data, size);

    //Create p-sample
    for (int i = 0; i < p; ++i) 
    {
        localPSample[i] = data[i * size / p];
        recDisplacements[i] = i * p;
        recCounts[i] = p;
    }

    *compTime += MPI::Wtime() - compTime2;
    commTime2 = MPI::Wtime();

    //Send all p-samples to proc 1
    MPI_Gatherv(localPSample, p, MPI_INT, allPSamples, recCounts, recDisplacements, MPI_INT, 0, MPI_COMM_WORLD);
    
    *commTime += MPI::Wtime() - commTime2;
    compTime2 = MPI::Wtime();

    //Sort all received samples and compute global p-sample
    if (id == 0)
    {
        heapsort(allPSamples, p * p);
        for (int i = 0; i < p; ++i)
        {
            localPSample[i] = allPSamples[i * p];
        }
    }
   
    *compTime += MPI::Wtime() - compTime2;
    commTime2 = MPI::Wtime();

    //Broadcast global p-sample
    MPI_Bcast(localPSample, p, MPI_INT, 0, MPI_COMM_WORLD);

    *commTime += MPI::Wtime() - commTime2;
    compTime2 = MPI::Wtime();

    //From this point on, localPSample contains the global p-sample 

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
    
    for (int i = 0; i < p-1; ++i)
    {
        bucketSize[i] = bucketLocation[i+1] - bucketLocation[i];
    }
    bucketSize[p-1] = size - bucketLocation[p-1];

    setupCountsAndDisplacements(sendCounts, recCounts, sendDisplacements, recDisplacements, 1, 1, -1, -1, p);

    *compTime += MPI::Wtime() - compTime2;
    commTime2 = MPI::Wtime();

    //Send bucket i to proc i
    MPI_Alltoallv(bucketSize, sendCounts, sendDisplacements, MPI_INT, allBuckets, recCounts, recDisplacements, MPI_INT, MPI_COMM_WORLD);

    *commTime += MPI::Wtime() - commTime2;
    compTime2 = MPI::Wtime();

    sortedDataSize = 1;
    for (int i = 0; i < p; ++i)
    {
        sendCounts[i] = bucketSize[i];
        recCounts[i] = allBuckets[i];
        sendDisplacements[i] = bucketLocation[i];
        recDisplacements[i] = sortedDataSize - 1; 
        sortedDataSize += allBuckets[i];
    }
    sortedDataSize--;

    *compTime += MPI::Wtime() - compTime2;
    commTime2 = MPI::Wtime();

    MPI_Alltoallv(data, sendCounts, sendDisplacements, MPI_INT, sortedData, recCounts, recDisplacements, MPI_INT, MPI_COMM_WORLD);

    *commTime += MPI::Wtime() - commTime2;
    compTime2 = MPI::Wtime();

    //Resort locally
    heapsort(sortedData, sortedDataSize);

    //Post processing: Array balancing
    //need n/p items per processor

    setupCountsAndDisplacements(sendCounts, recCounts, sendDisplacements, recDisplacements, 1, 1, 0, -1, p);
   
    *compTime += MPI::Wtime() - compTime2;
    commTime2 = MPI::Wtime();

    //Distribute current bucket size of each processor
    MPI_Alltoallv(&sortedDataSize, sendCounts, sendDisplacements, MPI_INT, allBuckets, recCounts, recDisplacements, MPI_INT, MPI_COMM_WORLD);

    *commTime += MPI::Wtime() - commTime2;
    compTime2 = MPI::Wtime();

    int cumulativeLeft = 0;
    int cumulativeRight = 0;
    int left = 0;
    int right = 0;

    //Calculate how much data to distribute to p-1 and p+1
    for (int i = 0; i < id; ++i)
    {
        cumulativeLeft += allBuckets[i];
    }
    cumulativeRight = cumulativeLeft + allBuckets[id] - 1;

    for (int i = 0; i < p; ++i)
    {
        left = i * size;
        right = (((i + 1) * size) - 1);
        balancingData[i] = 0;
        if ((cumulativeLeft <= left) && (left <= cumulativeRight) && (cumulativeRight <= right))
        {
            balancingData[i] = cumulativeRight - left + 1;
        }
        else if ((left <= cumulativeLeft) && (cumulativeLeft <= right) && (right <= cumulativeRight))
        {
            balancingData[i] = right - cumulativeLeft + 1;
        }
        else if ((left <= cumulativeLeft) && (cumulativeRight <= right))
        {
            balancingData[i] = cumulativeRight - cumulativeLeft + 1;
        }
        else if ((cumulativeLeft <= left) && (right <= cumulativeRight))
        {
            balancingData[i] = right - left + 1;
        }
    }

    setupCountsAndDisplacements(sendCounts, recCounts, sendDisplacements, recDisplacements, 1, 1, -1, -1, p);

    *compTime += MPI::Wtime() - compTime2;
    commTime2 = MPI::Wtime();

    MPI_Alltoallv(balancingData, sendCounts, sendDisplacements, MPI_INT, balancedData, recCounts, recDisplacements, MPI_INT, MPI_COMM_WORLD);

    *commTime += MPI::Wtime() - commTime2;
    compTime2 = MPI::Wtime();

    left = 0;
    right = 0;

    for (int i = 0; i < p; ++i)
    {
        sendCounts[i] = balancingData[i];
        recCounts[i] = balancedData[i];
        sendDisplacements[i] = left;
        recDisplacements[i] = right;
        left += sendCounts[i];
        right += recCounts[i];
    }

    *compTime += MPI::Wtime() - compTime2;
    commTime2 = MPI::Wtime();

    MPI_Alltoallv(sortedData, sendCounts, sendDisplacements, MPI_INT, finalData, recCounts, recDisplacements, MPI_INT, MPI_COMM_WORLD);

    *commTime += MPI::Wtime() - commTime2;
    compTime2 = MPI::Wtime();

    delete [] localPSample;
    delete [] allPSamples;
    delete [] allBuckets;
    delete [] recCounts;
    delete [] sendCounts;
    delete [] recDisplacements;
    delete [] sendDisplacements;
    delete [] bucketLocation;
    delete [] bucketSize;
    delete [] sortedData;
    delete [] balancingData;

    *compTime += MPI::Wtime() - compTime2;
}

int main(int argc, char *argv[])
{
    char inputFilePath[256];
    char outputFilePath[256];
    int n;
    int p;
    ifstream input;
    string buffer;
    int i;

    int size;

    int *data;
    int *finalData;

    int processors;
    int id;

    double startTime;
    double endTime;
    double computationTime = 0;
    double communicationTime = 0;

    MPI::Init(argc, argv); //  Initialize MPI.
    processors = MPI::COMM_WORLD.Get_size(); //  Get the number of processes.
    id = MPI::COMM_WORLD.Get_rank(); //  Get the individual process ID.

    sprintf(inputFilePath, "/tmp/fredericmarchand/input-%d.txt", id);
    sprintf(outputFilePath, "/tmp/fredericmarchand/output-%d.txt", id);
    
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
    finalData = new int[size + 1];

    while (getline(input, buffer))
    {
        data[i++] = atoi(buffer.c_str());
    }
    
    input.close();
    
    startTime = MPI::Wtime();
    deterministicSampleSort(data, size, finalData, p, id, &computationTime, &communicationTime);
    endTime = MPI::Wtime();

    cout << "Total Time for proc" << id << ": " << endTime - startTime << ", compTime: " << computationTime << ", commTime: " << communicationTime << endl;

    //Write to file
    ofstream outputFile (outputFilePath); 

    if (!outputFile.is_open())
    {
        cerr << "Could not open file " << outputFilePath << endl;
        return 1;   /// Error
    }

    for (int i = 0; i < size; ++i)
    {
        outputFile << finalData[i] << endl;
    }

    outputFile.close();

    // Terminate MPI.
    MPI::Finalize();

    delete [] data;
    delete [] finalData;

    return 0;
}
