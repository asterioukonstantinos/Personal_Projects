/***********************

Conway's Game of Life

Based on https://web.cs.dal.ca/~arc/teaching/CS4125/2014winter/Assignment2/Assignment2.html

************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <mpi.h>
#include <omp.h>

static int print_cells = 0;
static int print_world = 0;
static int root_process = 0;
static int random_world = 0;
static int root_rank = 0;
static int num_omp_threads = 0;

int world_size; // number of processes
int world_rank; // the rank of the process

static char *start_world[] = {
    /* Gosper glider gun */
    /* example from https://bitstorm.org/gameoflife/ */
    "..........................................",
    "..........................................",
    "..........................................",
    "..........................................",
    "..........................................",
    "..........................................",
    "........................OO.........OO.....",
    ".......................O.O.........OO.....",
    ".OO.......OO...........OO.................",
    ".OO......O.O..............................",
    ".........OO......OO.......................",
    ".................O.O......................",
    ".................O........................",
    "....................................OO....",
    "....................................O.O...",
    "....................................O.....",
    "..........................................",
    "..........................................",
    ".........................OOO..............",
    ".........................O................",
    "..........................O...............",
    "..........................................",
};

static void
world_init_fixed(int** current_world, int height, int width)
{
    int i, j;

    /* use predefined start_world */
    for (i = 1; i <= height; i++) {
        for (j = 1; j <= width; j++) {
            if ((i <= sizeof(start_world) / sizeof(char *)) &&
                (j <= strlen(start_world[i - 1]))) {
                current_world[i][j] = (start_world[i - 1][j - 1] != '.');
            } else {
                current_world[i][j] = 0;
            }
        }
    }
}

static void
world_init_random(int** current_world, int height, int width)
{
    int i, j;
    // Note that rand() implementation is platform dependent.
    // At least make it reprodible on this platform by means of srand()
    srand(1);

    for (i = 1; i <= height; i++) {
        for (j = 1; j <= width; j++) {
            float x = rand() / ((float)RAND_MAX + 1);
            if (x < 0.5) {
                current_world[i][j] = 0;
            } else {
                current_world[i][j] = 1;
            }
        }
    }
}

static void
world_print(int** current_world, int height, int width, int lb, int ub)
{
    int i, j, ub_last_proc;
     if(world_rank == 0){
        ub_last_proc += (ub - lb) + (height % world_size);
        for(int i = world_rank + 1; i < world_size; i++){
            if(i == world_size - 1)
                MPI_Recv(current_world[i*ub+1], (ub_last_proc-lb)*(width+2), MPI_INT, i, 0, MPI_COMM_WORLD, NULL);
            else            
                MPI_Recv(current_world[i*ub+1], (ub-lb)*(width+2), MPI_INT, i, 0, MPI_COMM_WORLD, NULL);
        }
    }
    else{
        MPI_Send(current_world[lb+1], (ub-lb)*(width+2), MPI_INT, 0, 0, MPI_COMM_WORLD);
        return;
    }

    for (i = 1; i <= height; i++) {
        for (j = 1; j <= width; j++) {
            if (current_world[i][j]) {
                printf("O");
            } else {
                printf(" ");
            }
        }
        printf("\n");
    }
}

static void
world_count(int** current_world, int width, int lb, int ub)
{
    int isum, total;
    int i, j;

    isum = 0;
    for (i = lb + 1; i <= ub; i++) {
        for (j = 1; j <= width; j++) {
            isum = isum + current_world[i][j];
        }
    }
    MPI_Reduce(&isum, &total, 1, MPI_INT, MPI_SUM, root_rank, MPI_COMM_WORLD);
    if(world_rank == 0)
        printf("Number of live cells = %d\n", total);
    return;
}

static int
world_count_total(int** current_world, int height, int width)
{
    int isum;
    int i, j;

    isum = 0;
    for (i = 1; i <= height; i++) {
        for (j = 1; j <= width; j++) {
            isum = isum + current_world[i][j];
        }
    }

    return isum;
}

/* Take world wrap-around into account: */
static void
world_border_wrap(int** current_world, int height, int width)
{
    int i, j;
    /* top-bottom boundary conditions */
    memcpy(current_world[0],current_world[height],sizeof(int)*(width+1));
    memcpy(current_world[height + 1],current_world[1],sizeof(int)*(width+1));
    
}

static void world_copy_cols(int** current_world, int lb, int ub, int width){
    int i, j;

    /* left-right boundary conditions */
    #pragma omp parallel for
    for (i = lb; i <= ub + 1; i++) {
        current_world[i][0] = current_world[i][width];
        current_world[i][width + 1] = current_world[i][1];
    }
}

static int world_cell_newstate(int** current_world, int row, int col)
{
    int row_m, row_p, col_m, col_p, nsum;
    int newval;

    // sum surrounding cells
    row_m = row - 1;
    row_p = row + 1;
    col_m = col - 1;
    col_p = col + 1;

    nsum = current_world[row_p][col_m] + current_world[row_p][col] + current_world[row_p][col_p]
         + current_world[row  ][col_m]                     + current_world[row  ][col_p]
         + current_world[row_m][col_m] + current_world[row_m][col] + current_world[row_m][col_p];

    switch (nsum) {
    case 3:
        // a new cell is born
        newval = 1;
        break;
    case 2:
        // nothing happens
        newval = current_world[row][col];
        break;
    default:
        // the cell, if any, dies
        newval = 0;
    }

    return newval;
}


// update board for next timestep
// height/width params are the base height/width
// excluding the surrounding 1-cell wraparound border
static void
world_timestep(int** current_world, int** next_world, int lb, int ub, int width)
{
    int row, col;
    
    // update board
    #pragma omp parallel for collapse(2)
    for (row = lb + 1; row <= ub; row++) {
        for (col = 1; col <= width; col++) {
                int row_m, row_p, col_m, col_p, nsum;
                int newval;

                // sum surrounding cells
                row_m = row - 1;
                row_p = row + 1;
                col_m = col - 1;
                col_p = col + 1;

                nsum = current_world[row_p][col_m] + current_world[row_p][col] + current_world[row_p][col_p]
                    + current_world[row  ][col_m]                     + current_world[row  ][col_p]
                    + current_world[row_m][col_m] + current_world[row_m][col] + current_world[row_m][col_p];


                //next_world[row][col] = ((nsum == 2) || (nsum == 3)) ? 1 : 0;
                switch (nsum) {
                case 3:
                    // a new cell is born
                    next_world[row][col] = 1;
                    break;
                case 2:
                    // nothing happens
                    next_world[row][col] = current_world[row][col];
                    break;
                default:
                    // the cell, if any, dies
                    next_world[row][col] = 0;
                }
        }
    }
}

static int **
alloc_2d_int_array(int nrows, int ncolumns)
{
    int **array;
    int i;

    /* version that keeps the 2d data contiguous, can help caching and slicing across dimensions */
    array = malloc(nrows * sizeof(int *));
    if (array == NULL) {
       fprintf(stderr, "out of memory\n");
       exit(1);
    }

    array[0] = malloc(nrows * ncolumns * sizeof(int));
    if (array[0] == NULL) {
       fprintf(stderr, "out of memory\n");
       exit(1);
    }

    for (i = 1; i < nrows; i++) {
	array[i] = array[0] + i * ncolumns;
    }

    return array;
}

static double
time_secs(void)
{
    struct timeval tv;

    if (gettimeofday(&tv, 0) != 0) {
        fprintf(stderr, "could not do timing\n");
        exit(1);
    }

    return tv.tv_sec + (tv.tv_usec / 1000000.0);
}

int
main(int argc, char *argv[])
{
    int n, nsteps, i;
    double start_time, end_time, elapsed_time;
    int bwidth, bheight;
    MPI_Request req;
    
    MPI_Init(NULL, NULL);      // initialize MPI environment
    // get the number of procs
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // get the Process-ID
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    /* Get Parameters */
    if (argc < 6) {
        fprintf(stderr, "Usage: %s width height steps print_world print_cells random_world num_omp_threads\n", argv[0]);
        exit(1);
    }
    bwidth = atoi(argv[1]);
    bheight = atoi(argv[2]);
    nsteps = atoi(argv[3]);
    print_world = atoi(argv[4]);
    print_cells = atoi(argv[5]);
    if(argc > 6){
        random_world = atoi(argv[6]);
        num_omp_threads = atoi(argv[7]);
    }

    omp_set_num_threads(num_omp_threads);

    int rows_per_proc =  bheight / world_size;
    int modulo = bheight % world_size;

    int lb = rows_per_proc*world_rank;
    int ub = lb + rows_per_proc;
    int ub_last_proc = 0;
    if(world_rank == (world_size - 1)){
        ub += modulo;
    }
    else if(world_rank == 0){
        ub_last_proc = (ub - lb) + modulo;
    }    
    /* initialize worlds, when allocating arrays, add 2 for ghost cells in both directorions */
    int** current_world = alloc_2d_int_array(bheight + 2, bwidth + 2);
    int** next_world = alloc_2d_int_array(bheight + 2, bwidth + 2);
    

    /*  initialize board */
    if (random_world) {
        world_init_random(current_world, bheight, bwidth);
    } else {
        world_init_fixed(current_world, bheight, bwidth);
    }

    if (print_world > 0) {
        printf("\ninitial world:\n\n");
        world_print(current_world, bheight, bwidth, lb, ub);
    }
    world_border_wrap(current_world, bheight, bwidth);
    int lower_proc = (world_rank == 0) ? (world_size - 1) : (world_rank - 1);
    int upper_proc = (world_rank == (world_size - 1)) ? 0 : (world_rank + 1);

    start_time = time_secs();
    int err = 0;
    /*  time steps */
    for (n = 0; n < nsteps; n++) {
        int** tmp_world;

        world_copy_cols(current_world, lb, ub, bwidth);
        world_timestep(current_world, next_world, lb, ub, bwidth);

        MPI_Isend(&next_world[lb+1][1], bwidth, MPI_INT, lower_proc, 0, MPI_COMM_WORLD ,&req);
        MPI_Isend(&next_world[ub][1], bwidth, MPI_INT, upper_proc, 1, MPI_COMM_WORLD, &req);

        MPI_Barrier(MPI_COMM_WORLD);
        
        MPI_Irecv(&next_world[lb][1], bwidth, MPI_INT, lower_proc, 1, MPI_COMM_WORLD, &req);
        MPI_Irecv(&next_world[ub+1][1], bwidth, MPI_INT, upper_proc, 0, MPI_COMM_WORLD, &req);
        
        // swap old and new worlds
        tmp_world = current_world;
        current_world = next_world;
        next_world = tmp_world;

        if (print_cells > 0 && (n % print_cells) == (print_cells - 1)) {
            world_count(current_world, bwidth, lb, ub);        }

        if (print_world > 0 && (n % print_world) == (print_world - 1)) {
            printf("\nafter time step %d:\n\n", n);
            world_print(current_world, bheight, bwidth, lb, ub);        
        }

    }

    end_time = time_secs();
    elapsed_time = end_time - start_time;
   
    /*  Iterations are done; sum the number of live cells */
    world_count(current_world, bwidth, lb, ub);
    if(world_rank == 0)
        fprintf(stderr, "Game of Life took %10.3f seconds\n", elapsed_time);

    free(current_world);
    free(next_world);
    MPI_Finalize();
    return 0;
}
