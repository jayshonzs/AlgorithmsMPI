#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>

int add(char *file1, char *file2, char *targetfile, int m, int n, int add)
{
	int rank, size;
	MPI_File fh1, fh2, fh3;
	int res1, res2;
	float data1[1024], data2[1024];
	float added[1024];
	int numrows = 0;
	int i = 0;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	numrows = m/size;
	if(rank < m%size)
	{
		numrows++;
	}
	
	MPI_File_open(MPI_COMM_WORLD, file1, MPI_MODE_RDWR, MPI_INFO_NULL, &fh1);
		
	res1 = loadmatrix_rows(&fh1, data1, numrows, rank, size, m, n);
	
	MPI_File_close(&fh1);
	
	MPI_File_open(MPI_COMM_WORLD, file2, MPI_MODE_RDWR, MPI_INFO_NULL, &fh2);
	
	res2 = loadmatrix_rows(&fh2, data2, numrows, rank, size, m, n);
	
	MPI_File_close(&fh2);
	
	//assert(res1==res2)
	
	if(add)
	{
		for(i = 0; i < res1; i++)
		{
			added[i] = data1[i] + data2[i];
			//printf("Proc %d: %f\n", rank, added[i]);
		}
	}
	else
	{
		for(i = 0; i < res1; i++)
		{
			added[i] = data1[i] - data2[i];
		}
	}
	
	MPI_File_open(MPI_COMM_WORLD, targetfile, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh3);
	
	savematrix_rows(&fh3, added, numrows, rank, size, m, n);
	
	MPI_File_close(&fh3);
	
	return 0;
}
