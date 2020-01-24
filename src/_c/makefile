omp:
	@gcc -fopenmp -O3 lu_openmp.c commons.c -o omp
pthread:
	@gcc -fopenmp lu_pthread.c commons.c -lm -lpthread -o pthread
clean:
	rm pthread omp
	clear