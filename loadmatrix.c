#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>

//load m*n matrix
int loadmatrix_rows(MPI_File *fh, float *rbuf, int numrows, int rank, int numtasks, int m, int n)
{
	MPI_Offset offset = 0;
	int i = 0, j = 0;
	MPI_Status status;
	MPI_Datatype rowtype;
	int result = 0;
	
	MPI_Type_contiguous(n, MPI_FLOAT, &rowtype);
	MPI_Type_commit(&rowtype);
	
	if(rank < m%numtasks)
	{
		offset = numrows*rank;
	}
	else if(rank == m%numtasks)
	{
		offset = (numrows+1)*rank;
	}
	else
	{
		offset = numrows*rank+m%numtasks;
	}
		
	MPI_File_set_view(*fh, 0, rowtype, rowtype, "native", MPI_INFO_NULL);
	
	result = MPI_File_read_at(*fh, offset, rbuf, numrows, rowtype, &status);
	if(result != MPI_SUCCESS)
	{
		printf("Proc %d read at %d error!\n", rank, offset);
	}
	
	/*
	for(i = 0; i < numrows; i++)
	{
		printf("Proc %d row %d: ", rank, i);
		for(j = 0; j < n; j++)
		{
			printf("%f, ", rbuf[i*n+j]);
		}
		printf("\n");
	}
	*/
	
	//MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Type_free(&rowtype);
	
	return n*numrows;
}

int loadmatrix_cols(MPI_File *fh, float *rbuf, int rank, int numtasks, int m, int n)
{
	float data[1024];
	MPI_Datatype darray;
	MPI_Status status;
	int gsizes[2] = {m, n};
	int distribs[2] = {MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_BLOCK};
	int dargs[2] = {MPI_DISTRIBUTE_DFLT_DARG, MPI_DISTRIBUTE_DFLT_DARG};
	int psizes[2] = {1, 4};
	int i, j;
	int count = 0;
	int cols = 0;
	
	cols = n/numtasks;
	if(rank < n%numtasks)
	{
		cols++;
	}
	
	MPI_Type_create_darray(numtasks, rank, 2, gsizes, distribs, dargs, psizes, MPI_ORDER_C, MPI_FLOAT, &darray);
	MPI_Type_commit(&darray);
	
	MPI_File_set_view(*fh, 0, MPI_FLOAT, darray, "native", MPI_INFO_NULL);
	
	MPI_File_read_all(*fh, data, m*cols, MPI_FLOAT, &status);
	
	MPI_Get_count(&status, MPI_FLOAT, &count);
	//cols = count/m;
	
	for(i = 0; i < cols; i++)
	{
		for(j = 0; j < m; j++)
		{
			rbuf[i*m+j] = data[j*cols+i];
		}
	}
	
	/*
	if(rank == 1)
	{
		for(i = 0; i < cols; i++)
		{
			printf("Proc %d row %d: ", rank, i);
			for(j = 0; j < m; j++)
			{
				printf("%f, ", rbuf[j*cols+i]);
			}
			printf("\n");
		}
	}
	*/
	
	MPI_Type_free(&darray);

	return cols;
}

int savematrix_rows(MPI_File *fh, float *data, int numrows, int rank, int numtasks, int m, int n)
{
	MPI_Datatype darray;
	MPI_Status status;
	int gsizes[2] = {m, n};
	int distribs[2] = {MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_BLOCK};
	int dargs[2] = {MPI_DISTRIBUTE_DFLT_DARG, MPI_DISTRIBUTE_DFLT_DARG};
	int psizes[2] = {4, 1};
	
	MPI_Type_create_darray(numtasks, rank, 2, gsizes, distribs, dargs, psizes, MPI_ORDER_C, MPI_FLOAT, &darray);
	MPI_Type_commit(&darray);
	
	MPI_File_set_view(*fh, 0, MPI_FLOAT, darray, "native", MPI_INFO_NULL);
	
	MPI_File_write_all(*fh, data, numrows*n, MPI_FLOAT, &status);
	
	MPI_Type_free(&darray);
	
	return 0;
}

int loadmatrix_whole(MPI_File *fh, float *data, int m, int n, int col)
{
	MPI_Status status;
	float **dataT;
	int i = 0, j = 0;

	MPI_File_set_view(*fh, 0, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL);
	
	MPI_File_read_all(*fh, data, m*n, MPI_FLOAT, &status);
	
	if(col)
	{
		dataT = (float**)malloc(n*sizeof(float*));
		for(i = 0; i < n; i++)
		{
			dataT[i] = (float*)malloc(m*sizeof(float));
		}
		
		for(i = 0; i < m; i++)
		{
			for(j = 0; j < n; j++)
			{
				dataT[j][i] = data[i*n+j];
			}
		}
		
		for(i = 0; i < n; i++)
		{
			for(j = 0; j < m; j++)
			{
				data[i*m+j] = dataT[i][j];
			}
		}
	}

	return 0;
}

int loadmatrix_cross_rows(MPI_File *fh, float *rbuf, int numrows, int rank, int numtasks, int m, int n)
{
	MPI_Offset offset = 0;
	int i = 0, j = 0;
	MPI_Status status;
	MPI_Datatype rowtype;
	MPI_Datatype filetype;
	int result = 0;
	
	MPI_Type_contiguous(n, MPI_FLOAT, &rowtype);
	MPI_Type_commit(&rowtype);
	
	MPI_Type_vector(numrows, 1, numtasks, rowtype, &filetype);
	MPI_Type_commit(&filetype);
	
	offset = rank;
	
	MPI_File_set_view(*fh, offset*n*sizeof(float), rowtype, filetype, "native", MPI_INFO_NULL);
	
	result = MPI_File_read_at(*fh, 0, rbuf, numrows, rowtype, &status);
	if(result != MPI_SUCCESS)
	{
		printf("Proc %d read at %d error!\n", rank, offset);
	}
	
	if(rank == 2)
	{
		for(i = 0; i < numrows; i++)
		{
			printf("Proc %d read row %d: ", rank, i);
			for(j = 0; j < n; j++)
			{
				printf("%f, ", rbuf[i*n+j]);
			}
			printf("\n");
		}
	}
	
	MPI_Type_free(&rowtype);
	MPI_Type_free(&filetype);
	
	return n*numrows;
}

int savematrix_cross_rows(MPI_File *fh, float *data, int numrows, int rank, int numtasks, int m, int n)
{
	MPI_Offset offset = 0;
	int i = 0, j = 0;
	MPI_Status status;
	MPI_Datatype rowtype;
	MPI_Datatype filetype;
	int result = 0;
	
	MPI_Type_contiguous(n, MPI_FLOAT, &rowtype);
	MPI_Type_commit(&rowtype);
	
	MPI_Type_vector(numrows, 1, numtasks, rowtype, &filetype);
	MPI_Type_commit(&filetype);
	
	offset = rank;
	
	MPI_File_set_view(*fh, offset*n*sizeof(float), rowtype, filetype, "native", MPI_INFO_NULL);
	
	result = MPI_File_write_at(*fh, 0, data, numrows, rowtype, &status);
	if(result != MPI_SUCCESS)
	{
		printf("Proc %d write at %d error!\n", rank, offset);
	}
	
	#if 1
	if(rank == 0)
	{
		for(i = 0; i < numrows; i++)
		{
			printf("Proc %d write row %d: ", rank, i);
			for(j = 0; j < n; j++)
			{
				printf("%f, ", data[i*n+j]);
			}
			printf("\n");
		}
	}
	#endif
	
	MPI_Type_free(&rowtype);
	MPI_Type_free(&filetype);
	
	return n*numrows;
}
