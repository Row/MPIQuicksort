########################################################################
# Makefile for MPIQuicksort Project
#
########################################################################

CC         =  mpicc
CCFLAGS    =  -O3
LIBS       =  -lmpi

qsort: qsort.c
	$(CC) $(CCFLAGS) -o qsort.out qsort.c $(LIBS)

nodebug: qsort.c
	$(CC) $(CCFLAGS) -DNDEBUG -o qsort.out qsort.c $(LIBS)