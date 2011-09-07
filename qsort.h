#ifndef MPI_QSORT_H_GUARD
#define MPI_QSORT_H_GUARD
#define ROOT_RESULT 666 
#define GROUP_EXCHANGE 667
    void print_array(int *a, int length);
    void print_sort_error(int *a, int length);
    void array_fill_rand(int *a, int length);
    void array_fill(int *a, int length);
    void testSuit();
    int is_sorted(int *a, int length);
    int is_range(int *a, int length);
    int is_equal(int a[], int b[], int length);
    int divide(int A[], int length, int pivot);
    void merge(int m, int n, int A[], int B[], int C[]);
    int mpi_qsort(int *IN[], int length, MPI_Comm procs_comm);
    int partition(int A[], int length);
    int compare(const void * a, const void * b);
#endif
