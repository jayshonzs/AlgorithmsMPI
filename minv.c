#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>

int loadmatrix_cross_rows(MPI_File *fh, float *rbuf, int numrows, int rank, int numtasks, int m, int n);
int savematrix_cross_rows(MPI_File *fh, float *data, int numrows, int rank, int numtasks, int m, int n);

int inv(char *srcfile, char *invfile, int n)
{
	float rbuf[1024];
	float f[512];
	//float rf[512];
	int rank, size, numrows;
	MPI_File fh;
	int i, j, k, w;
	int row = 0;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	numrows = n/size;
	if(rank < n%size)
	{
		numrows++;
	}
	
	MPI_File_open(MPI_COMM_WORLD, srcfile, MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
	
	loadmatrix_cross_rows(&fh, rbuf, numrows, rank, size, 8, 8);
	
	MPI_File_close(&fh);
	
	for(i = 0; i < numrows; i++)
	{
		for(j = 0; j < size; j++)
		{
			if(rank == j)
			{
				row = i*size + j;
				rbuf[i*n+row] = 1.0/rbuf[i*n+row];
				
				for(k = 0; k < n; k++)
				{
					if(k != row)
					{
						rbuf[i*n+k] = rbuf[i*n+k]*rbuf[i*n+row];
					}
				}
				
				for(k = 0; k < n; k++)
				{
					f[k] = rbuf[i*n+k];
				}
				
				//scatter main row
				//wrong:MPI_Scatter(f, n, MPI_FLOAT, rf, n, MPI_FLOAT, rank, MPI_COMM_WORLD);
				MPI_Bcast(f, n, MPI_FLOAT, rank/*root process*/, MPI_COMM_WORLD);
				
				//printf("Proc %d dest %d f: %f, %f, %f, %f, %f, %f, %f, %f\n",
				//	rank, j, f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7]);
			}
			else
			{
				row = i*size + j;
				//recv main row
				//wrong:MPI_Scatter(f, n, MPI_FLOAT, rf, n, MPI_FLOAT, j, MPI_COMM_WORLD);
				MPI_Bcast(f, n, MPI_FLOAT, j/*root process*/, MPI_COMM_WORLD);
				
				//printf("Proc %d dest %d f: %f, %f, %f, %f, %f, %f, %f, %f\n",
				//	rank, j, f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7]);
			}
			
			//MPI_Barrier(MPI_COMM_WORLD);
			//exit(-1);
			
			if(rank != j)
			{
				for(k = 0; k < numrows; k++)
				{
					for(w = 0; w < n; w++)
					{
						if(w != row)
						{
							rbuf[k*n+w] = rbuf[k*n+w] - f[w]*rbuf[k*n+row];
						}
					}
				}
				
				for(k = 0; k < numrows; k++)
				{
					rbuf[k*n+row] = -f[row]*rbuf[k*n+row];
				}
			}
			else
			{
				for(k = 0; k < numrows; k++)
				{
					if(k != i)
					{
						for(w = 0; w < n; w++)
						{
							if(w != row)
							{
								rbuf[k*n+w] = rbuf[k*n+w] - f[w]*rbuf[k*n+row];
							}
						}
					}
				}
				
				for(k = 0; k < numrows; k++)
				{
					if(k != i)
					{
						rbuf[k*n+row] = -f[row]*rbuf[k*n+row];
					}
				}
			}
		}
	}
	
	MPI_File_open(MPI_COMM_WORLD, invfile, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	savematrix_cross_rows(&fh, rbuf, numrows, rank, size, n, n);
	MPI_File_close(&fh);
	
	return 0;
}
