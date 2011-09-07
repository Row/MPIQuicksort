/**********************************************************************
 * Quicksort in MPI/C with 2^n processors.
 *
 * 2011 Jon Borglund
 *
 * mpirun -np 8 qsort.out
 * mpirun -hostfile nodes -mca plm_rsh_agent rsh -np 64 qsort.out
 **********************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "qsort.h"

void print_array(int *a, int length)
{
    int i;
    for (i = 0; i < length; i++)
        printf("%d ", a[i]);

    printf(", len: %d\n", length);
}

void print_sort_error(int *a, int length)
{
    int i;
    for (i = 1; i < length; i++)
        if (a[i - 1] > a[i])
            printf("error: a[%d] = %d, a[%d] = %d\n", i-1, a[i-1], i, a[i]); ;

}

int is_sorted(int *a, int length)
{
    int i;
    for (i = 1; i < length; i++)
        if (a[i - 1] > a[i])
            return 0;

    return 1;
}

int is_range(int *a, int length)
{
    int i;
    for (i = 0; i < length; i++)
        if (a[i] != i+1)
            return 0;

    return 1;
}

int is_equal(int a[], int b[], int length)
{
    int i;
    for (i = 0; i < length; i++)
        if (a[i] != b[i])
            return 0;

    return 1;
}

void array_fill_rand(int *a, int length)
{
    int i;
    for (i = 0; i < length; i++)
        a[i] = rand();
}

void array_fill(int *a, int length)
{
    int i;
    for (i = 0; i < length; i++)
        a[i] = length - i;
}

void testSuit()
{
    int A[7] = {1, 2, 3, 4, 5, 6};
    int B[7] = {1, 2, 2, 3, 4, 5, 6};
    int C[7] = {2, 2, 3, 1, 4, 5, 6};
    int D[7] = {1, 2, 3, 4, 5, 6};
    int E[7] = {0, 1, 2, 3, 4, 5, 6};
    int F[7] = {0, 1, 2, 2, 3, 4, 5};
    int G[7] = {0, 2, 2, 3, 4, 5, 6};
    int H[7] = {0, 4, 8, 16};
    assert(is_range(A, 6));
    assert(!is_range(B, 7));

    assert(is_sorted(A, 6));
    assert(is_sorted(B, 7));
    assert(!is_sorted(C, 7));

    assert(is_equal(A, D, 6));
    assert(!is_equal(B, D, 6));

    assert(is_equal(A, D, 6));

    assert(divide(E, 7, 3) == 4);
    assert(divide(E, 7, -1) == 0);
    assert(divide(E, 7, 0) == 1);
    assert(divide(E, 7, 10) == 7);
    assert(divide(F, 7, 2) == 4);
    assert(divide(F, 7, 1) == 2);
    assert(divide(G, 7, 1) == 1);
    assert(divide(G, 7, 2) == 3);
    assert(divide(H, 4, 0) == 1);
    assert(divide(H, 4, 1) == 1);
    assert(divide(H, 4, 8) == 3);
    assert(divide(H, 4, 10) == 3);
    assert(divide(H, 4, 5) == 2);
}

/****************************************************************************
 **    MPI Quicksort
 ****************************************************************************/
int partition(int A[], int length)
{
    return (A[0] + A[length - 1]) / 2;
}

int divide(int A[], int length, int pivot)
{
    int divider = length/2;
    while (divider < length && A[divider] <= pivot)
        divider++;
    while (divider > 0 && A[divider - 1] > pivot)
        divider--;

    return divider;
}

void merge(int m, int n, int A[], int B[], int C[])
{
    int i, j, k, p;
    i = j = k = 0;
    while (i < m && j < n) {
        if (A[i] <= B[j]) {
            C[k] = A[i];
            i++;
        } else {
            C[k] = B[j];
            j++;
        }
        k++;
    }
    if (i < m) {
        for (p = i; p < m; p++) {
            C[k] = A[p];
            k++;
        }
    } else {
        for (p = j; p < n; p++) {
            C[k] = B[p];
            k++;
        }
    }
}

int mpi_qsort(int *IN[], int length, MPI_Comm procs_comm)
{
    int number_of_processors, myrank, color, pivot, divider;
    int dest, offset, chunksize, recv_length, newlength;
    int *A, *B, *C;

    MPI_Comm group_comm;
    MPI_Status status;
    MPI_Request request;
    MPI_Comm_size(procs_comm, &number_of_processors);
    MPI_Comm_rank(procs_comm, &myrank);

    // Return length of new array
    if (number_of_processors < 2)
        return length;

    A = *IN;

    // Select pivot and distr
    pivot = partition(A, length);
    MPI_Bcast(&pivot, 1, MPI_INT, 0, procs_comm);

    // Set divide index
    divider = divide(A, length, pivot);

    // Divide into sets
    color = myrank / (number_of_processors / 2);

    if (color) {
        // Right side send left
        dest = myrank - number_of_processors / 2;
        offset = 0;
        chunksize = divider;
    } else {
        // Left side Send right
        dest = myrank + number_of_processors / 2;
        offset = divider;
        chunksize = length - divider;
    }

    // Send
    MPI_Isend(&A[offset], chunksize, MPI_INT, dest, GROUP_EXCHANGE, procs_comm, &request);

    // Get length of message and receive
    MPI_Probe(dest, GROUP_EXCHANGE, procs_comm, &status);
    MPI_Get_count(&status, MPI_INT, &recv_length);
    B = malloc(recv_length*sizeof(*B));
    MPI_Recv(B, recv_length, MPI_INT, dest, GROUP_EXCHANGE, procs_comm, &status);

    // Merge old part and newly received part
    newlength = recv_length + length - chunksize;
    offset = offset ? 0 : divider;
    C = malloc(newlength*sizeof(*C));
    merge(length - chunksize, recv_length, &A[offset], B, C);

    // Wait for Isend to finish
    MPI_Wait(&request, &status);

    // Deallocate and set IN pointer to new merged array
    free(A);
    free(B);
    *IN = C;

    // Split communicator and do recursion
    MPI_Comm_split(procs_comm, color, myrank, &group_comm);
    return mpi_qsort(IN, newlength, group_comm);
}

int compare(const void * a, const void * b)
{
  return (*(int*)a - *(int*)b);
}

int main(int argc, char *argv[])
{
    int n = argc > 1 ? atoi(argv[1]) : 100;

    double time;

    //Number of processes
    int number_of_processors;

    //Process rank/id
    int myrank;

    MPI_Status status;
    MPI_Comm procs;
    MPI_Request request;

    /* Initialize MPI               */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processors);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    //Size of each local processors array
    int local_size = n/number_of_processors;

    // Buffers
    int *local_array, *root_array;
    local_array = malloc(local_size*sizeof(*local_array));

    if (myrank == 0) {
        printf("MPI Quicksort: %d processes, %d elements.\n", number_of_processors, n);
        root_array = malloc(n*sizeof(*local_array));
        array_fill_rand(root_array, n);
        time = MPI_Wtime();
    }

    // Distribute the root array
    MPI_Scatter(root_array, local_size, MPI_INT, local_array, local_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Sort locally
    qsort(local_array, local_size, sizeof(*local_array), compare);

    // Peform parallel sorting
    local_size = mpi_qsort(&local_array, local_size, MPI_COMM_WORLD);

    assert(is_sorted(local_array, local_size));

    // Send result to root process
    MPI_Isend(local_array, local_size, MPI_INT, 0, ROOT_RESULT, MPI_COMM_WORLD, &request);

    //Gather result
    if (myrank == 0) {
        int recv_length, offset, p;
        offset = 0;
        for(p = 0; p < number_of_processors; p++) {
            MPI_Recv(&root_array[offset], n, MPI_INT, p, ROOT_RESULT, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_INT, &recv_length);
            offset += recv_length;
        }

        time = MPI_Wtime() - time;

        printf("%f\n", time, n);

        // Some tests
        //print_array(root_array, n);
        assert(is_sorted(root_array, n));
        //print_sort_error(root_array, n);
        testSuit();
        //assert(is_range(root_array, n));

    }

    // Wait for send to root process
    MPI_Wait(&request, &status);

    MPI_Finalize();
    return 0;
}
