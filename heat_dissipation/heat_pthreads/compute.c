#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "compute.h"
#include <pthread.h>
//#define IMAGE


void calc_weights(double* conductivity, double* diagonalweight, double* neighborweight, size_t rows, size_t cols){
	double diagonal = 1 /(1 + sqrt(2));
	double neighbor = sqrt(2) / (1 + sqrt(2));
	for(int i = 1; i < (rows-1); i++){
		for(int j = 0; j < cols; j++){
            diagonalweight[i*cols + j] = ((1 - conductivity[(i-1)*cols + j])*diagonal)/4;
            neighborweight[i*cols + j] = ((1 - conductivity[(i-1)*cols + j])*neighbor)/4;
		}
	}
}


void* do_compute(void* a)
{   
    extern struct parameters p;
    extern struct results r;
    extern double *told;
    extern double *tnext;
    extern pthread_barrier_t  barrier;
    struct timespec before, after;
	size_t rows = p.N + 2, cols = p.M, iter, row, col, temp_pos, cond_pos, left, right, bottom, top;
    double temp, diff, maxdiff, tmin, tmax, tsum;

    int* tid_p = a;   // cast  arg to int  pointer
    int tid = *((int*)a);     
    int rowsperthread = p.N/p.nthreads;
	double *tempold = told, *tempnew = tnext, *temptmp;
	double *neighborweight = calloc(rows * cols, sizeof(double));
	double *diagonalweight = calloc(rows * cols, sizeof(double));
    calc_weights((double *)p.conductivity,diagonalweight,neighborweight,rows,cols);	// calculate weights of each element
    int lower_bound = tid*rowsperthread + 1;
    int upper_bound = lower_bound + rowsperthread;
    if((p.N % p.nthreads != 0) && (tid == (p.nthreads - 1))){
        upper_bound += (p.N % p.nthreads);
    }
    
    clock_gettime(CLOCK_MONOTONIC, &before);
 	for(iter = 0; iter < p.maxiter; iter++){
        maxdiff = 0;
		for(row = lower_bound; row < upper_bound; row++){		// for each row
			for(col = 0; col < cols; col++){	// for each column
				temp_pos = row*cols + col;
                cond_pos = temp_pos - cols; // cond_pos = (row-1)*cols + col; // the position in the matrix. Replace the whole expression-more readable

                left = ((col == 0) ? (cols - 1) : (col - 1));
                right = ((col + 1 == cols) ? 0 : (col + 1));

                bottom = row + 1;
                top = row - 1;
				tempnew[temp_pos] = p.conductivity[cond_pos] * tempold[temp_pos] +
				        diagonalweight[temp_pos] * (tempold[top*cols + left] + tempold[top*cols + right] + tempold[bottom*cols + left] + tempold[bottom*cols + right]) +
				        neighborweight[temp_pos] * (tempold[row*cols + left] + tempold[row*cols + right] + tempold[top*cols + col] + tempold[bottom*cols + col]);

                diff = fabs(tempnew[temp_pos] - tempold[temp_pos]);
                if (diff > maxdiff) maxdiff = diff;

			}
		}
        temptmp = tempold;
        tempold = tempnew;
        tempnew = temptmp;    
        pthread_barrier_wait(&barrier);
		if (maxdiff < p.threshold) break;
 	}
    if(tid != 0){
        return a;
    } 
    temptmp = tempold;
    tempold = tempnew;
    tempnew = temptmp;

    tmin = DBL_MAX;
    tmax = DBL_MIN;
    tsum = 0;
    for(row = 1; row < (rows-1); row++) {
        for (col = 0; col < cols; col++) {
            temp = tempnew[row * cols + col];
            tsum += temp;
            if (temp > tmax) tmax = temp;
            if (temp < tmin) tmin = temp;
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &after);
    r.niter = iter;
 	r.tmin = tmin;
    r.tmax = tmax;
 	r.maxdiff = maxdiff;
    r.tavg = tsum / (double) (p.N * p.M);
    r.time = (double)(after.tv_sec - before.tv_sec) +
              (double)(after.tv_nsec - before.tv_nsec) / 1e9;
	free(neighborweight);
	free(diagonalweight);
	
    
    return a;
}
