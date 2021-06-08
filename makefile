CC= gcc
CFLAGS=-Wall -fopenmp
LMATH=-lm
EXCLUDE=-ffunction-sections
OFLAG =-Ofast -march=native -funroll-loops -ftree-vectorize
DEBUGG = #-fopt-info-vec -g -pg

eigen: eigenvalues.c
	$(CC) $(CFLAGS) $(OFLAG) $(DEBUGG) -o eigen eigenvalues.c $(LMATH) $(EXCLUDE)

clean:
	rm -f eigen 
