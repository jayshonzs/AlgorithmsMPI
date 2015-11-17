#include <stdio.h>
#include <mpi.h>
#include <memory.h>
#include <stdlib.h>
#include <unistd.h>

int loadmatrix_rows(MPI_File *fh, float *rbuf, int numrows, int rank, int numtasks, int m, int n);
int loadmatrix_cols(MPI_File *fh, float *rbuf, int rank, int numtasks, int m, int n);

float dot(float *a, float *b, int length, int rank)
{
	float ret = 0;
	int i = 0;
	
	if(rank == -1)
	{
		printf("dot length=%d\n", length);
	}
	
	for(i = 0; i < length; i++)
	{
		ret += a[i]*b[i];
		if(rank == -1)
		{
			printf("%f, ", a[i]);
		}
	}
	if(rank == -1)
	{
		printf("\n\n");
	}
	
	return ret;
}

int multiply(char *leftfile, char *rightfile, char *targetfile, int left_m, int left_n, int right_m, int right_n)
{
	MPI_Status status;
	MPI_File fh1, fh2, fh3;
	int rank, size;
	float leftdata[1024*4], rightdata[1024*4];
	float targetdata[1024*4];
	int res = 0, cols = 0, col_offset = 0;
	int rcols = 0;
	float rrightdata[1024*4];
	int numrows = 0;
	int i = 0, j = 0, k = 0;
	int left = 0, right = 0;
	
	//assert(left_n==right_m);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	numrows = left_m/size;
	if(rank < left_m%size)
	{
		numrows++;
	}
	
	left = rank - 1;
	right = rank + 1;
	if(left < 0)
	{
		left = size-1;
	}
	if(right >= size)
	{
		right = 0;
	}
	
	printf("Proc %d numrows: %d\n", rank, numrows);
	//MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_File_open(MPI_COMM_WORLD, leftfile, MPI_MODE_RDWR, MPI_INFO_NULL, &fh1);
	res = loadmatrix_rows(&fh1, leftdata, numrows, rank, size, left_m, left_n);
	MPI_File_close(&fh1);
	
	/*
	if(rank == 0)
	{
		for(i = 0; i < numrows; i++)
		{
			printf("Proc %d row %d: ", rank, i);
			for(j = 0; j < left_n; j++)
			{
				printf("%f, ", leftdata[i*left_n+j]);
			}
			printf("\n");
		}
	}
	*/
	
	MPI_File_open(MPI_COMM_WORLD, rightfile, MPI_MODE_RDWR, MPI_INFO_NULL, &fh2);
	if(right_n >= size)
	{
		cols = loadmatrix_cols(&fh2, rightdata, rank, size, right_m, right_n);
		
		MPI_Scan(&cols, &col_offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		col_offset -= cols;
		printf("Proc %d cols=%d, col_offset = %d\n", rank, cols, col_offset);
	}
	else
	{
		cols = right_n;
		loadmatrix_whole(&fh2, rightdata, right_m, right_n, 1);
	}
	MPI_File_close(&fh2);
	
	/*
	if(rank == 1)
	{
		for(i = 0; i < cols; i++)
		{
			printf("Proc %d col %d: ", rank, i);
			for(j = 0; j < right_m; j++)
			{
				printf("%f, ", rightdata[i*right_m+j]);
			}
			printf("\n");
		}
	}
	*/
	
	printf("Proc %d res:%d, cols:%d\n", rank, res, cols);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(right_n >= size)
	{
		for(i = 0; i < size; i++)
		{
			
			for(j = 0; j < numrows; j++)
			{
				for(k = 0; k < cols; k++)
				{
					targetdata[j*right_n+col_offset+k] = dot(&leftdata[j*left_n], &rightdata[k*right_m], left_n, rank);
					if(rank == -1)
					{
						printf("Proc %d: j=%d, k=%d, %f\n", rank, j, k, targetdata[j*right_n+col_offset+k]);
					}
				}
			}
			
			col_offset += cols;
			if(col_offset == right_n)
			{
				col_offset = 0;
			}
			
			MPI_Sendrecv(&cols,
							1,
							MPI_INT,
							left,
							123,
							&rcols,
							1,
							MPI_INT,
							right,
							123,
							MPI_COMM_WORLD,
							&status);
			
			MPI_Sendrecv(rightdata,
							cols*right_m,
							MPI_FLOAT,
							left,
							124,
							rrightdata,
							rcols*right_m,
							MPI_FLOAT,
							right,
							124,
							MPI_COMM_WORLD,
							&status);
			
			cols = rcols;
			memcpy(rightdata, rrightdata, 1024*4*sizeof(float));
		}
	}
	else
	{
		for(j = 0; j < numrows; j++)
		{
			for(k = 0; k < cols; k++)
			{
				targetdata[j*right_n+k] = dot(&leftdata[j*left_n], &rightdata[k*right_m], left_n, rank);
				printf("Proc %d : %f\n", rank, targetdata[j*right_n+k]);
			}
		}
	}
	
	if(rank == -1)
	{
		for(i = 0; i < numrows; i++)
		{
			printf("Proc %d row %d: ", rank, i);
			for(j = 0; j < right_n; j++)
			{
				printf("%f, ", targetdata[i*right_n+j]);
			}
			printf("\n");
		}
	}
	
	MPI_File_open(MPI_COMM_WORLD, targetfile, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh3);
	cols = savematrix_rows(&fh3, targetdata, numrows, rank, size, left_m, right_n);
	MPI_File_close(&fh3);

	return 0;
}
