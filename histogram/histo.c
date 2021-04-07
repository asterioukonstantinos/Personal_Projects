#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include <pthread.h>

struct info_image
{	
	int* image;
	int rows;
	int columns;
	int threads_num;
};

struct info_image* img;
int** histo;

void die(const char *msg){
    if (errno != 0) 
        perror(msg);
    else
        fprintf(stderr, "error: %s\n", msg);
    exit(1);
}   

void generate_image(int num_rows, int num_cols, int * image){
    for (int i = 0; i < num_cols * num_rows; ++i)
    {
        image[i] = 1; //255 + 1 for num bins
    }
}

void read_image(const char * image_path, int num_rows, int num_cols, int * image){
	char format[3];
    FILE *f;
    unsigned imgw, imgh, maxv, v;
    size_t i;

	printf("Reading PGM data from %s...\n", image_path);

	if (!(f = fopen(image_path, "r"))) die("fopen");

	fscanf(f, "%2s", format);
    if (format[0] != 'P' || format[1] != '2') die("only ASCII PGM input is supported");
    
    if (fscanf(f, "%u", &imgw) != 1 ||
        fscanf(f, "%u", &imgh) != 1 ||
        fscanf(f, "%u", &maxv) != 1) die("invalid input");

    if (imgw != num_cols || imgh != num_rows) {
        fprintf(stderr, "input data size (%ux%u) does not match cylinder size (%zux%zu)\n",
                imgw, imgh, num_cols, num_rows);
        die("invalid input");
    }

    for (i = 0; i < num_cols * num_rows; ++i)
    {
        if (fscanf(f, "%u", &v) != 1) die("invalid data");
        image[i] = ((int)v * 255) / maxv; //255 for num bins
    }
    fclose(f);
}

void print_histo(int * histo){
	for (int i = 0; i < 256; ++i)
	{	
		if(i != 0 && (i % 10 == 0)) {
            printf("\n");
        }
		printf("%d ", histo[i]);
	}
    printf("\n");
}

void print_image(int num_rows, int num_cols, int * image){
	int index = 0;
	for (int i = 0; i < num_rows; ++i){	
		for (int j = 0; j < num_cols; ++j){
	        index = i * num_cols + j;
			printf("%d ", image[index]);
		}
	}
	printf("\n");
}

void* histogram(void* id){
	int index = 0;
    int thread_id = *((int*)id);
    int rows_per_thread = img->rows/img->threads_num;
	int lower_bound = *((int*)id)*rows_per_thread;
	int upper_bound = lower_bound + rows_per_thread;
    if((img->rows%img->threads_num != 0) && (thread_id == (img->threads_num - 1))){
        upper_bound += img->rows%img->threads_num;
    }
    for (int i = lower_bound; i < upper_bound; ++i){	
		for (int j = 0; j < img->columns; ++j){
	        index = i * img->columns + j;
			histo[thread_id][img->image[index]]++;
		}
	}
	return id;
}

int main(int argc, char *argv[]){
    int c;
    int seed = 42;
    const char *image_path = 0;
    image_path ="../../../../images/pat1_100x150.pgm";
    int gen_image = 0;
    int debug = 0;
 	img = (struct info_image*)malloc(sizeof(struct info_image));
    int num_rows = 150;
    int num_cols = 100;
    int num_threads = 1;
    struct timespec before, after;
    void* result;
    int* num;
    int* final = (int *) calloc(256, sizeof(int));

    /* Read command-line options. */
    while((c = getopt(argc, argv, "s:i:rp:n:m:g")) != -1) {
        switch(c) {
            case 's':
                seed = atoi(optarg);
                break;
            case 'i':
            	image_path = optarg;
            	break;
            case 'r':
            	gen_image = 1;
            	break;
            case 'p':
                num_threads = atoi(optarg);
                break;
            case 'n':
            	num_rows = strtol(optarg, 0, 10);
            	break;
            case 'm':
				num_cols = strtol(optarg, 0, 10);
				break;
			case 'g':
				debug = 1;
				break;
            case '?':
                fprintf(stderr, "Unknown option character '\\x%x'.\n", optopt);
                return -1;
            default:
                return -1;
        }
    }
    img->rows = num_rows;
    img->columns = num_cols;
    img->threads_num = num_threads; 
    img->image = (int *) malloc(sizeof(int) * num_cols * num_rows);
    pthread_t  thread_ids[img->threads_num];
    //printf("I am here\n");

    /* Seed such that we can always reproduce the same random vector */
    if (gen_image){
    	srand(seed);
    	generate_image(img->rows, img->columns, img->image);
    }else{
    	read_image(image_path,img->rows, img->columns, img->image);
    }

    histo = (int**)malloc(sizeof(int*)*img->threads_num);
    for(int i = 0; i < img->threads_num; i++){
    	histo[i] = (int*)calloc(256,sizeof(int));
    }


    clock_gettime(CLOCK_MONOTONIC, &before);
    /* Do your thing here */
    for (int i=0; i<img->threads_num; i++) {
    	num = (int*)  malloc( sizeof( int));
    	*num = i;
    	pthread_create( &thread_ids[i], NULL ,&histogram , num);
    	//printf("Created all thread[%d]\n",i);
    }

    /* Do your thing here */

    for (int i=0; i<img->threads_num; i++) {
    	pthread_join( thread_ids[i], &result );
    	free( result );
    }

    for(int i = 0; i < 256; i++){
    	for(int j = 0; j < img->threads_num; j++){
    		final[i] += histo[j][i];
    	}
    }

    clock_gettime(CLOCK_MONOTONIC, &after);
    double time = (double)(after.tv_sec - before.tv_sec) +
                  (double)(after.tv_nsec - before.tv_nsec) / 1e9;

    if (debug){
    	print_histo(final);
    }
    
    for(int i = 0; i < img->threads_num; i++){
    	free(histo[i]);
    }
    free(histo);
    printf("Histo took: % .6e seconds \n", time);

}
