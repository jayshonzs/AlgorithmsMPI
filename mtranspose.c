#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>

int loadmatrix_cols(MPI_File *fh, float *rbuf, int rank, int numtasks, int m, int n);

/*
 * m*n matrix transposing
 */
int transpose(char *srcfile, char *targetfile, int m, int n)
{
	float rbuf[1024];
	//float data[1024];
	int size, rank;
	MPI_File fh;
	int cols = 0, i, j;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	MPI_File_open(MPI_COMM_WORLD, srcfile, MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
	
	cols = loadmatrix_cols(&fh, rbuf, rank, size, 67, 8);
	
	MPI_File_close(&fh);
	
	MPI_File_open(MPI_COMM_WORLD, targetfile, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	
	savematrix_rows(&fh, rbuf, cols, rank, size, 8, 67);
	
	MPI_File_close(&fh);
		
	return 0;
}
