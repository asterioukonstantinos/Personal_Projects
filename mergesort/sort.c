#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include "omp.h"
#include "sort.h"


void CopyArray(int* A,long iBegin,long iEnd,int* B)
{
    for(long k = iBegin; k < iEnd; k++)
        B[k] = A[k];
}

//  Left source half is A[ iBegin:iMiddle-1].
// Right source half is A[iMiddle:iEnd-1   ].
// Result is            B[ iBegin:iEnd-1   ].
void TopDownMerge(int* A,long iBegin,long iMiddle,long iEnd,int* B)
{
    long i = iBegin, j = iMiddle;
    long size = iEnd - iBegin;
    // While there are elements in the left or right runs...
    #pragma omp simd
    for (long k = iBegin; k < iEnd; k++) {
        // If left run head exists and is <= existing right run head.
        if (i < iMiddle && (j >= iEnd || A[i] <= A[j])) {
            B[k] = A[i];
            i = i + 1;
        } else {
            B[k] = A[j];
            j = j + 1;
        }
    }
    
    // printf("\n");
    // for(long k = iBegin; k < iEnd; k++){
    //     printf("%d  ",B[k]);
    // }
    // printf( "\nThe Thread %d of %d finished. iBegin = %ld, iEnd = %ld, iMiddle = %ld\n\n",omp_get_thread_num(),omp_get_num_threads(), iBegin, iEnd, iMiddle);
}

// Sort the given run of array A[] using array B[] as a source.
// iBegin is inclusive; iEnd is exclusive (A[iEnd] is not in the set).
void TopDownSplitMerge(int* B,long iBegin,long iEnd,int* A)
{
    if(iEnd - iBegin <= 1)                       // if run size == 1
        return;                                 //   consider it sorted
    // split the run longer than 1 item into halves
    long iMiddle = (iEnd + iBegin) / 2;              // iMiddle = mid point
    
    #pragma omp taskgroup
    {
        #pragma omp task firstprivate (iBegin, iMiddle, iEnd) mergeable untied if(iMiddle - iBegin > 1000)\
                         shared(A,B)
        {
            TopDownSplitMerge(A, iBegin,  iMiddle, B);  // sort the left  run
            #pragma omp taskyield
            //printf( "Thread %d of %d finished.\n",omp_get_thread_num (),omp_get_num_threads ());
            // merge the resulting runs from array B[] into A[]
        }
        #pragma omp task firstprivate (iBegin, iMiddle, iEnd) mergeable untied if(iEnd - iMiddle > 1000)\
                         shared(A,B)
        {
            TopDownSplitMerge(A, iMiddle, iEnd, B);  // sort the right run
            #pragma omp taskyield
            //printf( "Thread %d of %d finished.\n",omp_get_thread_num (),omp_get_num_threads ());
            // merge the resulting runs from array B[] into A[]
        }
    }       
   
    // merge the resulting runs from array B[] into A[]
    
    TopDownMerge(B, iBegin, iMiddle, iEnd, A);
}

// Array A[] has the items to sort; array B[] is a work array.
void mergesort(int* A, long n)
{
    int* B = (int*)malloc(sizeof(int)*n);
    CopyArray(A, 0, n, B);           // one time copy of A[] to B[]
    #pragma omp parallel
    {
        #pragma omp single nowait
        TopDownSplitMerge(B, 0, n, A);   // sort data from B[] into A[]
    }
    free(B);
}






