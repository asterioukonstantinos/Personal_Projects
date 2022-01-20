# When gol-par.c has been created, uncomment the next line:
# all: gol-seq gol-par

gol-seq: gol-seq.c
	gcc -Wall -O3 -o gol-seq gol-seq.c -lm

gol-par: gol-par.c
	mpicc -Wall -O3 -mavx2 -fopenmp -o gol-par gol-par.c -lm

gol-par-bonus: gol-par-bonus.c 
	mpicc -Wall -O3 -mavx2 -fopenmp -o gol-par-bonus gol-par-bonus.c -lm

clean:
	rm -f *.o gol-seq gol-par *~ *core
