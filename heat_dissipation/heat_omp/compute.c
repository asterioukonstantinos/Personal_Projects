#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "compute.h"
#include "omp.h"


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

void initialize(double* told,double* tinit, size_t rows, size_t cols){

	for(int i = 0; i < cols; i++){
		told[i] = tinit[i];
	}
	for(int i = 0; i < (rows-2); i++){
		for(int j = 0; j < cols; j++){
			told[(i+1)*cols + j] = tinit[i*cols + j];
		}
	}
	for(int i = 0; i < cols; i++){
		told[(rows-1)*cols + i] = tinit[(rows-3)*cols + i];
	}
}

void do_compute(const struct parameters* p, struct results *r)
{
    struct timespec before, after;
	size_t rows = p->N + 2, cols = p->M, iter, row, col, temp_pos, cond_pos, left, right, bottom, top;
    double *conductivity = calloc((p->N*p->M),sizeof(double));
    for(int i = 0; i < (p->N*p->M); i++){
        conductivity[i] = p->conductivity[i];
    }
	double *told = calloc(rows * cols, sizeof(double));	// allocate matrix for olds values
	double *tnext = calloc(rows * cols, sizeof(double)); // allocate matrix for new values
	double *tempold = told, *tempnew = tnext, *temptmp;
	double *neighborweight = calloc(rows * cols, sizeof(double));
	double *diagonalweight = calloc(rows * cols, sizeof(double));
	double temp, diff, maxdiff, tmin, tmax, tsum;

    initialize(told,(double *)p->tinit,rows,cols);	// initiliaze allocated array
    initialize(tnext,(double *)p->tinit,rows,cols);	// initiliaze allocated array
    calc_weights((double *)conductivity,diagonalweight,neighborweight,rows,cols);	// calculate weights of each element
    
    printf("Number of cores -> %d\n",omp_get_num_procs());
	printf("COMPUTATION START\n");
    
    int block = ((p->N*p->M) / omp_get_num_procs());
    tmin = DBL_MAX;
    tmax = DBL_MIN;
    clock_gettime(CLOCK_MONOTONIC, &before);
 	for(iter = 0; iter < p->maxiter; iter++){	// for that many iterations
        maxdiff = 0;
        tsum = 0;
        
        #pragma  omp  parallel  for  schedule(static,block) collapse(2)\
                            reduction( max:maxdiff)\
                            shared( tempnew, tempold)\
                            private(temp_pos, cond_pos, left, right, bottom, top)
		for(row = 1; row < (rows-1); row++){		// for each row
			for(col = 0; col < cols; col++){	// for each column
				temp_pos = row*cols + col;
                cond_pos = temp_pos - cols; // the position in the matrix. Replace the whole expression-more readable

                left = ((col == 0) ? (cols - 1) : (col - 1));
                right = ((col + 1 == cols) ? 0 : (col + 1));
                bottom = row + 1;
                top = row - 1;
				tempnew[temp_pos] =  diagonalweight[temp_pos] * (tempold[top*cols + left] + tempold[top*cols + right] + tempold[bottom*cols + left] + tempold[bottom*cols + right]) 
                + neighborweight[temp_pos] * (tempold[row*cols + left] + tempold[row*cols + right] + tempold[top*cols + col] + tempold[bottom*cols + col])
                + conductivity[cond_pos] * tempold[temp_pos];
                if ((fabs(tempnew[temp_pos] - tempold[temp_pos])) > maxdiff) 
                    maxdiff = fabs(tempnew[temp_pos] - tempold[temp_pos]);
			}
        }
        temptmp = tempold;
		tempold = tempnew;
		tempnew = temptmp;
        if(maxdiff < p->threshold){
            break;
        }
    }

     // ---- DO NOT REMOVE / FIXES ACCURACY ----------
    temptmp = tempold;
    tempold = tempnew;
    tempnew = temptmp;
    
    
    tmin = DBL_MAX;
    tmax = DBL_MIN;
    #pragma  omp  parallel  for  schedule(static) num_threads(2)\
                              reduction( +: tsum)\
                              reduction( max:tmax)\
                              reduction( min:tmin)\
                              shared(tempnew)
    for(row = 1; row < (rows-1); row++) {
        for (col = 0; col < cols; col++) {
            tsum += tempnew[temp_pos];
            if (tempnew[row * cols + col] > tmax) tmax = tempnew[row * cols + col];
            if (tempnew[row * cols + col] < tmin) tmin = tempnew[row * cols + col];
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &after);
    printf("COMPUTATION FINISH\n");
    r->niter = iter;
 	r->tmin = tmin;
    r->tmax = tmax;
 	r->maxdiff = maxdiff;
    r->tavg = tsum / (double) (rows * cols);
    r->time = (double)(after.tv_sec - before.tv_sec) +
              (double)(after.tv_nsec - before.tv_nsec) / 1e9;
    free(conductivity);
	free(neighborweight);
	free(diagonalweight);
	free(told);
	free(tnext);
}
