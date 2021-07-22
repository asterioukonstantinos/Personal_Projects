#include "compute.h"
#define _GNU_SOURCE
#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>

struct parameters p;
struct results r;
double *told;
double *tnext;
pthread_barrier_t   barrier; // the barrier synchronization object


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


int main(int argc, char **argv)
{
	read_parameters(&p, argc, argv);
	size_t rows = p.N + 2, cols = p.M;
	told = calloc(rows * cols, sizeof(double));	// allocate matrix for olds values
	tnext = calloc(rows * cols, sizeof(double)); // allocate matrix for new values
	initialize(told,(double *)p.tinit,rows,cols);	// initiliaze allocated array
    initialize(tnext,(double *)p.tinit,rows,cols);	// initiliaze allocated array
    pthread_t  thread_ids[p.nthreads];
    pthread_barrier_init (&barrier, NULL, p.nthreads);
	int i, *num;
	void *result;
	for (i=0; i<p.nthreads; i++) 
	{
		num = (int*)malloc(sizeof(int));
		*num = i;
		pthread_create( &thread_ids[i], NULL ,&do_compute, num);
	}
	for (i=0; i<p.nthreads; i++) 
	{
		pthread_join( thread_ids[i], &result );
		free(result);
	}

    report_results(&p, &r);
    free(told);
	free(tnext);
    return 0;
}

