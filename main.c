#include <stdio.h>
#include <mpi.h>

int kmeans(int K, double **centers);
int loadmatrix_rows(MPI_File *fh, int rank, int numtasks, int m, int n);
int loadmatrix_cols(MPI_File *fh, float *rbuf, int rank, int numtasks, int m, int n);
int add(char *file1, char *file2, char *targetfile, int m, int n, int add);
int multiply(char *leftfile, char *rightfile, char *targetfile, int left_m, int left_n, int right_m, int right_n);
int loadmatrix_cross_rows(MPI_File *fh, float *rbuf, int numrows, int rank, int numtasks, int m, int n);
int inv(char *srcfile, char *invfile, int n);

int ridged_regression(char *xfile, char *yfile, int m, int n, float lambda);

int main(int argc, char **argv)
{
	double centers[3][2];
	float rbuf[1024];
	
	MPI_Init(0, 0);
	
	//kmeans(3, (double**)centers);
	
	//int numrows = 0;
	//int size, rank;
	//MPI_File fh;
	//MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//MPI_Comm_size(MPI_COMM_WORLD, &size);
	//MPI_File_open(MPI_COMM_WORLD, "binary_X.dat", MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
	
	//test load rows
	/*
	loadmatrix_rows(&fh, rank, size, 67, 8);
	*/
	
	//test load columns
	//loadmatrix_cols(&fh, rbuf, rank, size, 67, 8);
	
	//MPI_File_close(&fh);
	
	//transpose("binary_X.dat", "binary_XT.dat", 67, 8);
	
	//add("binary_X.dat", "binary_X.dat", "binary_X_ADD.dat", 67, 8, 1);
	//add("binary_X.dat", "binary_X.dat", "binary_X_SUB.dat", 67, 8, 0);
	
	//multiply("binary_XT.dat", "binary_X.dat", "binary_XT_dot_X.dat", 8, 67, 67, 8);
	//multiply("binary_XT.dat", "binary_Y.dat", "binary_XT_dot_Y.dat", 8, 67, 67, 1);
	
	//test load cross rows
	//numrows = 67/size;
	//if(rank < 67%size)
	//{
	//	numrows++;
	//}
	//loadmatrix_cross_rows(&fh, rbuf, numrows, rank, size, 67, 8);
	//MPI_File_close(&fh);
	
	//test inv
	//inv("binary_XT_dot_X.dat", "binary_inv.dat", 8);
	
	//ridged test
	ridged_regression("binary_X.dat", "binary_Y.dat", 67, 8, 23);
	
	MPI_Finalize();
	
	return 0;
}
